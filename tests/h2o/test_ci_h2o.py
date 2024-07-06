
import numpy as np
from fcipy.driver import Driver

def test_ci_h2o():
    driver = Driver.from_gamess(logfile="h2o-Re.log", onebody="onebody.inp", twobody="twobody.inp", nfrozen=0)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants(file="psi_det_1000")

    # run CI calculation
    driver.run_ci(nroot=5)

    # check the results
    assert np.allclose(driver.total_energy[0], -76.22471425)

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
    # cisolver = fci.FCI(mf)
    # print('E(FCI) = %.12f' % cisolver.kernel(frozen=1)[0])

    driver = Driver.from_pyscf(mf, nfrozen=0)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants(method="fci", target_irrep="A1")

    # run CI calculation
    driver.run_ci(nroot=1)

    # check the results
    assert np.allclose(driver.total_energy[0], -75.012009000900)

if __name__ == "__main__":
    #test_ci_h2o()
    test_fci_h2o()
