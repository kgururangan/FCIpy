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

if __name__ == "__main__":
    test_ci_h2o()
