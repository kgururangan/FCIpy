"""Closed-shell CCSDTQ calculation for the symmetrically stretched
H2O molecule with R(OH) = 2Re, where Re = 1.84345 bohr, described
using the Dunning DZ basis set.
Reference: Mol. Phys, 115, 2860 (2017)."""

import numpy as np
from fcipy.driver import Driver

def test_ci_h2o():
    driver = Driver.from_gamess(logfile="h2o-Re.log", onebody="onebody.inp", twobody="twobody.inp", nfrozen=0)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants("psi_det_1000")

    # run CI calculation
    driver.run_ci(nroot=5)

if __name__ == "__main__":
    test_ci_h2o()
