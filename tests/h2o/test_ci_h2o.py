import numpy as np
from fcipy.driver import Driver

from pathlib import Path
TEST_DATA_DIR = str(Path(__file__).parents[0].absolute() / "data")

def test_ci_h2o():

    # Test that FCI solver can reproduce the results from QP2 calculations using
    # general lists of determinants obtained via the CIPSI algorithm
    n_det_max = ["100", "1000", "10000"]
    expected_energy = [-76.131590287, -76.224714253, -76.23571330]

    for ndet, e_expected in zip(n_det_max, expected_energy):

        driver = Driver.from_gamess(logfile=TEST_DATA_DIR + "/h2o-Re.log", 
                                    onebody=TEST_DATA_DIR + "/onebody.inp", 
                                    twobody=TEST_DATA_DIR + "/twobody.inp", 
                                    nfrozen=0)
        driver.system.print_info()

        # read the determinants
        driver.load_determinants(file=TEST_DATA_DIR + f"/psi_det_{ndet}")

        # run CI calculation
        driver.run_ci(nroot=1, opt=True, herm=False)

        # check the results
        assert np.allclose(driver.total_energy[0], e_expected)

    _, state_energy_compute = driver.compute_energy(herm=False)
    
    for e1, e2 in zip(state_energy_compute, driver.total_energy):
        assert np.allclose(e1, e2, rtol=1.0e-08, atol=1.0e-08)

if __name__ == "__main__":
    test_ci_h2o()
