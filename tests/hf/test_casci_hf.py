import numpy as np
from pyscf import gto, scf, mcscf
from fcipy.driver import Driver

if __name__ == "__main__":

    mol = gto.M(atom='''H 0 0 0; F 0 0 1.5''', 
                spin=0, 
                basis='6-31g', 
                unit='angstrom', 
                verbose=0, 
                symmetry="C2v"
    )
    mf = scf.RHF(mol).run()

    # Total orbital space
    nels = mol.nelectron
    norb = mf.mo_coeff.shape[1]

    # Set up active space: CAS(6,4)
    nelcas_a = 3
    nelcas_b = 3
    nact = 4

    nelcas = nelcas_a + nelcas_b
    ncore = (nels - nelcas)//2
    nvirt = norb - nact - ncore

    cas = (nelcas, nact)

    #
    # create an FCI solver based on the SCF object
    #
    driver = Driver.from_pyscf(mf, nfrozen=2, ndelete=5)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants(max_excit_rank=-1, target_irrep="A1")

    # run CI calculation in the active space
    driver.run_ci(nroot=1, prtol=0.01, herm=True)

    # Run CASCI
    mc = mcscf.CASCI(mf, cas[1], cas[0])
    mc.run()

    assert np.allclose(driver.total_energy[0], mc.e_tot, atol=1.0e-08, rtol=1.0e-08)
