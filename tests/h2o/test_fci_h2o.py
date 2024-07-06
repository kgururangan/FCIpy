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

    driver = Driver.from_pyscf(mf, nfrozen=0)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants(method="fci", target_irrep="A1")

    # run CI calculation
    driver.run_ci(nroot=1)

    # check the results
    assert np.allclose(driver.total_energy[0], cisolver.kernel(frozen=0)[0])

if __name__ == "__main__":
    test_fci_h2o()