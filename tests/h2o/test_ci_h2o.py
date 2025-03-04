import numpy as np
from fcipy.driver import Driver

def test_ci_h2o():

    # Test that FCI solver can reproduce the results from QP2 calculations using
    # general lists of determinants obtained via the CIPSI algorithm
    n_det_max = ["100", "1000", "10000"]
    expected_energy = [-76.131590287, -76.224714253, -76.23571330]

    for ndet, e_expected in zip(n_det_max, expected_energy):

        driver = Driver.from_gamess(logfile="h2o-Re.log", onebody="onebody.inp", twobody="twobody.inp", nfrozen=0)
        driver.system.print_info()

        # read the determinants
        driver.load_determinants(file=f"psi_det_{ndet}")

        # run CI calculation
        driver.run_ci(nroot=1)

        # check the results
        assert np.allclose(driver.total_energy[0], e_expected)

if __name__ == "__main__":
    test_ci_h2o()
