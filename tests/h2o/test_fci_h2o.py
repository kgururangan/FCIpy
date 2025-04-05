import numpy as np
from fcipy.driver import Driver

def test_fci_h2o():
    from pyscf import gto, scf, fci

    geom = [["H", (0, 1.515263, -1.058898)],
            ["H", (0, -1.515263, -1.058898)],
            ["O", (0.0, 0.0, -0.0090)]]

    mol = gto.M(atom=geom, basis="sto-3g", spin=0, charge=0, unit="Bohr", symmetry="C2V")
    mf = scf.RHF(mol)
    mf.kernel()

    #
    # create an FCI solver based on the SCF object
    #
    cisolver = fci.FCI(mf)

    driver = Driver.from_pyscf(mf, nfrozen=1)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants(max_excit_rank=-1, target_irrep="A1")

    # run CI calculation
    driver.run_ci(nroot=1, prtol=0.01)
    #driver.diagonalize_hamiltonian()
    print(driver.total_energy[0])

    # compute the 1-RDM
    #driver.one_e_density_matrix()

    #print("Natural occupation numbers:")
    #for i, n in enumerate(driver.nat_occ_num[0]):
    #    print(f"Orbital {i + 1}: {n}")

    # check the results
    assert np.allclose(driver.total_energy[0], cisolver.kernel(frozen=1)[0])

if __name__ == "__main__":
    test_fci_h2o()
