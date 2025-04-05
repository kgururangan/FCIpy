import numpy as np
from pyscf import ao2mo, symm
from fcipy.system import System

eVtohartree = 0.036749308136649
hartreetoeV = 1.0/(eVtohartree)

def load_pyscf_integrals(meanfield, nfrozen=0, ndelete=0):
    """Builds the System and Integral objects using the information contained within a PySCF
    mean-field object for a molecular system.

    Arguments:
    ----------
    meanFieldObj : Object -> PySCF SCF/mean-field object
    nfrozen : int -> number of frozen electrons
    Returns:
    ----------
    system: System object
    integrals: Integral object
    """
    molecule = meanfield.mol
    nelectrons = molecule.nelectron
    mo_coeff = meanfield.mo_coeff
    norbitals = mo_coeff.shape[1]
    nuclear_repulsion = molecule.energy_nuc()

    system = System(
        nelectrons,
        norbitals,
        molecule.spin + 1,  # PySCF mol.spin returns 2S, not S
        nfrozen,
        ndelete=ndelete,
        point_group=molecule.symmetry,
        orbital_symmetries = [x.upper() for x in symm.label_orb_symm(molecule, molecule.irrep_name, molecule.symm_orb, mo_coeff)],
        charge=molecule.charge,
        nuclear_repulsion=nuclear_repulsion,
        mo_energies=meanfield.mo_energy,
        mo_occupation=meanfield.mo_occ,
    )

    # Perform AO-to-MO transformation (using mf.get_hcore() allows this to work with scalar 1e X2C models, for instance)
    e1int = np.einsum(
        "pi,pq,qj->ij", mo_coeff, meanfield.get_hcore(), mo_coeff, optimize=True
    )
    e1int = np.asfortranarray(e1int)

    # two-body integrals
    e2int = np.transpose(
        np.reshape(ao2mo.kernel(molecule, mo_coeff, compact=False), 4 * (norbitals,)),
        (0, 2, 1, 3)
    )
    e2int = np.asfortranarray(e2int)
    # Check that the HF energy calculated using the integrals matches the PySCF result
    hf_energy = calc_hf_energy(e1int, e2int, system)
    hf_frozen = calc_hf_frozen_energy(e1int, e2int, system)
    #frozen_energy = hf_energy - hf_unfrozen
    hf_energy += nuclear_repulsion

    if not np.allclose(hf_energy, meanfield.energy_tot(), atol=1.0e-06, rtol=0.0):
        raise RuntimeError("Integrals don't match mean field energy")

    system.reference_energy = hf_energy
    system.frozen_energy = hf_frozen
    if system.nfrozen > 0:
        K = np.einsum("ipiq->pq", e2int[:nfrozen, nfrozen:, :nfrozen, nfrozen:])
        J = np.einsum("ipqi->pq", e2int[:nfrozen, nfrozen:, nfrozen:, :nfrozen])
        e1int[nfrozen:, nfrozen:] += 2*K - J

    return system, e1int[nfrozen:, nfrozen:], e2int[nfrozen:, nfrozen:, nfrozen:, nfrozen:]

def load_gamess_integrals(
    logfile=None,
    fcidump=None,
    onebody=None,
    twobody=None,
    nfrozen=0,
    ndelete=0,
    multiplicity=None,
    data_type=np.float64,
):
    from cclib.io import ccread

    data = ccread(logfile)
    # unable to change nelectrons attribute
    nelectrons = data.nelectrons
    # if multiplicity is provided, use that; otherwise, read from GAMESS logfile
    if multiplicity is None:
        multiplicity = data.mult
    # Read nelectrons & norbitals directly from fcidump in order to accomodate ECP GAMESS runs
    if fcidump:
        nmo_fcidump, nelectron_fcidump, ms2_fcidump = load_system_params_from_fcidump(fcidump)
        data.nmo = nmo_fcidump
        if data.nelectrons != nelectron_fcidump:
            print("\n   =============================================================================================")
            print("      WARNING: NELEC in GAMESS logfile and FCIDUMP do not match! Using NELEC listed in FCIDUMP.")
            print("      (Beware, GAMESS FCIDUMP prints wrong NELEC and MS2 for open-shell systems!)")
            print("   =============================================================================================\n")
        nelectrons = nelectron_fcidump

    system = System(
        nelectrons,
        data.nmo,
        multiplicity,
        nfrozen,
        ndelete=ndelete,
        point_group=get_point_group(logfile),
        orbital_symmetries=[x.upper() for x in data.mosyms[0]],
        charge=data.charge,
        nuclear_repulsion=get_nuclear_repulsion(logfile),
        mo_energies=[x * eVtohartree for x in data.moenergies[0]],
    )

    # Load using onebody and twobody direct integral files
    if fcidump is None and onebody is not None and twobody is not None:
        e1int = load_onebody_integrals(onebody, system, data_type)
        nuclear_repulsion, e2int = load_twobody_integrals(twobody, system, data_type)
    # Load from FCIDUMP file
    elif fcidump is not None:
        e1int, e2int, nuclear_repulsion = load_integrals_from_fcidump(fcidump, system)

    assert np.allclose(
        nuclear_repulsion, system.nuclear_repulsion, atol=1.0e-06, rtol=0.0
    )
    system.nuclear_repulsion = nuclear_repulsion

    # Check that the HF energy calculated using the integrals matches the GAMESS result
    hf_energy = calc_hf_energy(e1int, e2int, system)
    hf_frozen = calc_hf_frozen_energy(e1int, e2int, system)
    #frozen_energy = hf_energy - hf_unfrozen
    hf_energy += system.nuclear_repulsion
    assert np.allclose(
        hf_energy, get_reference_energy(logfile), atol=1.0e-06, rtol=0.0
    )
    system.reference_energy = hf_energy
    system.frozen_energy = hf_frozen
    if system.nfrozen > 0:
        K = np.einsum("ipiq->pq", e2int[:nfrozen, nfrozen:, :nfrozen, nfrozen:])
        J = np.einsum("ipqi->pq", e2int[:nfrozen, nfrozen:, nfrozen:, :nfrozen])
        e1int[nfrozen:, nfrozen:] += 2*K - J

    return system, e1int[nfrozen:, nfrozen:], e2int[nfrozen:, nfrozen:, nfrozen:, nfrozen:]

def get_reference_energy(gamess_logfile):

    with open(gamess_logfile, "r") as f:
        for line in f.readlines():
            if all(s in line.split() for s in ["FINAL", "ROHF", "ENERGY", "IS"]) or all(
                s in line.split() for s in ["FINAL", "RHF", "ENERGY", "IS"]
            ):
                hf_energy = float(line.split()[4])
                break
    return hf_energy

def get_nuclear_repulsion(gamess_logfile):

    with open(gamess_logfile, "r") as f:
        for line in f.readlines():
            if all(
                s in line.split()
                for s in ["THE", "NUCLEAR", "REPULSION", "ENERGY", "IS"]
            ):
                e_nuclear = float(line.split()[-1])
                break
    return e_nuclear

def get_point_group(gamess_logfile):
    """Dumb way of getting the point group from GAMESS log files.

    Arguments:
    ----------
    gamessFile : str -> Path to GAMESS log file
    Returns:
    ----------
    point_group : str -> Molecular point group"""
    point_group = "C1"
    flag_found = False
    with open(gamess_logfile, "r") as f:
        for line in f.readlines():
            if flag_found and point_group != "CI" and point_group != "CS":
                order = line.split()[-1]
                if len(point_group) == 3:
                    point_group = point_group[0] + order + point_group[2]
                if len(point_group) == 2:
                    point_group = point_group[0] + order
                if len(point_group) == 1:
                    point_group = point_group[0] + order
                break
            if "THE POINT GROUP OF THE MOLECULE IS" in line:
                point_group = line.split()[-1]
                flag_found = True
    if point_group == 'C0':
        point_group = 'C1'
    return point_group

def load_onebody_integrals(onebody_file, system, data_type):
    """This function reads the onebody.inp file from GAMESS
    and returns a numpy matrix.

    Parameters
    ----------
    filename : str
        Path to onebody integral file
    sys : dict
        System information dict

    Returns
    -------
    e1int : ndarray(dtype=float, shape=(norb,norb))
        Onebody part of the bare Hamiltonian in the MO basis (Z)
    """
    norb = system.norbitals + system.nfrozen - system.ndelete
    e1int = np.zeros((norb, norb), dtype=data_type, order="F")
    try:
        with open(onebody_file) as f_in:
            lines = f_in.readlines()
            ct = 0
            for i in range(norb):
                for j in range(i + 1):
                    val = float(lines[ct].split()[0])
                    e1int[i, j] = val
                    e1int[j, i] = val
                    ct += 1
    except IOError:
        print("Error: {} does not appear to exist.".format(onebody_file))
    return e1int

def load_twobody_integrals(twobody_file, system, data_type):
    """This function reads the twobody.inp file from GAMESS
    and returns a numpy matrix.

    Parameters
    ----------
    filename : str
        Path to twobody integral file
    sys : dict
        System information dict

    Returns
    -------
    e_nn : float
        Nuclear repulsion energy (in hartree)
    e2int : ndarray(dtype=float, shape=(norb,norb,norb,norb))
        Twobody part of the bare Hamiltonian in the MO basis (V)
    """
    try:
        norb = system.norbitals + system.nfrozen - system.ndelete
        # initialize numpy array
        e2int = np.zeros((norb, norb, norb, norb), dtype=data_type, order="F")
        # open file
        with open(twobody_file) as f_in:
            # loop over lines
            for line in f_in:
                # split fields and parse
                fields = line.split()
                indices = tuple(map(int, fields[:4]))
                val = float(fields[4])
                # check whether value is nuclear repulsion
                # fill matrix otherwise
                if sum(indices) == 0:
                    e_nn = val
                else:
                    indices = tuple(i - 1 for i in indices)
                    e2int[indices] = val
        # convert e2int from chemist notation (ia|jb) to
        # physicist notation <ij|ab>
        e2int = np.einsum("iajb->ijab", e2int)
    except IOError:
        print("Error: {} does not appear to exist.".format(twobody_file))
    return e_nn, e2int

def load_system_params_from_fcidump(fcidump):
    """This function parses and returns the values printed in the top line of
    the FCIDUMP file, which lists the number of orbitals (NORB), the number of
    electrons (NELEC), and 2*S_z (MS2) for the system."""
    with open(fcidump, "r") as f:
        # Holds ['NORB=#', 'NELEC=#', 'MS2=#']
        firstline = f.readline().strip("\n").split(',')
        for entry in firstline:
            if "NORB" in entry:
                # Get norbitals
                norbitals = int(entry.split("=")[1])
            if "NELEC" in entry:
                # Get nelectrons
                nelectrons = int(entry.split("=")[1])
            if "MS2" in entry:
                # Get MS2 (note: this is not going to be correct for open shells in GAMESS)
                ms2 = int(entry.split("=")[1])
    return norbitals, nelectrons, ms2

def load_integrals_from_fcidump(fcidump, system):
    """This function reads the FCIDUMP file to obtain the onebody and twobody
    integrals as well as nuclear repulsion energy.

    Parameters
    ----------
    fcidump : str
        Path to FCIDUMP file
    system : System object
        System object

    Returns
    -------
    e1int : ndarray(dtype=float, shape=(norb,norb))
        Onebody part of the bare Hamiltonian in the MO basis (Z)
    e2int : ndarray(dtype=float, shape=(norb,norb,norb,norb))
        Twobody part of the bare Hamiltonian in the MO basis (V)
    e_nn : float
        Nuclear repulsion energy (in hartree)
    """
    norb = system.norbitals + system.nfrozen - system.ndelete
    e1int = np.zeros((norb, norb), order="F")
    e2int = np.zeros((norb, norb, norb, norb), order="F")

    with open(fcidump) as f:
        for ct, line in enumerate(f.readlines()):
            if ct < 4: continue
            L = line.split()
            Cf = float(L[0].replace("D", "E")) # GAMESS FCIDUMP uses old-school D instead of E for scientific notation
            p = int(L[1]) - 1
            q = int(L[3]) - 1
            r = int(L[2]) - 1
            s = int(L[4]) - 1
            if q != -1 and s != -1: # twobody term
                e2int[p, q, r, s] = Cf
                e2int[r, q, p, s] = Cf
                e2int[p, s, r, q] = Cf
                e2int[r, s, p, q] = Cf
                e2int[q, p, s, r] = Cf
                e2int[q, r, s, p] = Cf
                e2int[s, p, q, r] = Cf
                e2int[s, r, q, p] = Cf
            elif q == -1 and s == -1 and p != -1: # onebody term
                e1int[p, r] = Cf
                e1int[r, p] = Cf
            else: # nuclear repulsion
                e_nn = Cf

    return e1int, e2int, e_nn

def calc_hf_energy(e1int, e2int, system):

    occ_a = slice(0, system.noccupied_alpha + system.nfrozen)
    occ_b = slice(0, system.noccupied_beta + system.nfrozen)

    e1a = np.einsum("ii->", e1int[occ_a, occ_a])
    e1b = np.einsum("ii->", e1int[occ_b, occ_b])
    e2a = 0.5 * (
        np.einsum("ijij->", e2int[occ_a, occ_a, occ_a, occ_a])
        - np.einsum("ijji->", e2int[occ_a, occ_a, occ_a, occ_a])
    )
    e2b = np.einsum("ijij->", e2int[occ_a, occ_b, occ_a, occ_b])
    e2c = 0.5 * (
        np.einsum("ijij->", e2int[occ_b, occ_b, occ_b, occ_b])
        - np.einsum("ijji->", e2int[occ_b, occ_b, occ_b, occ_b])
    )

    hf_energy = e1a + e1b + e2a + e2b + e2c

    return hf_energy

def calc_hf_frozen_energy(e1int, e2int, system):

    occ_a = slice(0, system.nfrozen)
    occ_b = slice(0, system.nfrozen)
    Nf = system.nfrozen

    # e1a = np.einsum("ii->", e1int[occ_a, occ_a])
    # e1b = np.einsum("ii->", e1int[occ_b, occ_b])
    # e2a = 0.5 * (
    #     np.einsum("ijij->", e2int[occ_a, occ_a, occ_a, occ_a])
    #     - np.einsum("ijji->", e2int[occ_a, occ_a, occ_a, occ_a])
    # )
    # e2b = np.einsum("ijij->", e2int[occ_a, occ_b, occ_a, occ_b])
    # e2c = 0.5 * (
    #     np.einsum("ijij->", e2int[occ_b, occ_b, occ_b, occ_b])
    #     - np.einsum("ijji->", e2int[occ_b, occ_b, occ_b, occ_b])
    # )

    hf_frozen = 2 * np.trace(e1int[:Nf, :Nf]) + \
          2 * np.einsum('ijij', e2int[:Nf, :Nf, :Nf, :Nf]) - \
          np.einsum('ijji', e2int[:Nf, :Nf, :Nf, :Nf])

    #hf_frozen = e1a + e1b + e2a + e2b + e2c

    return hf_frozen
