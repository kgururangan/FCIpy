import numpy as np
from fcipy.lib import ci

def energy_ci_state(det, coef, e1int, e2int, noa, nob, num_alpha, num_beta, herm):

    e_diag = ci.ci.calc_diagonal(det, noa, nob, e1int, e2int)

    if herm:
        sigma = ci.ci.calc_sigma_opt(det, e_diag, coef, num_alpha, num_beta, e1int, e2int, noa, nob)

    else:
        sigma = ci.ci.calc_sigma_nonhermitian_opt(det, e_diag, coef, num_alpha, num_beta, e1int, e2int, noa, nob)

    energy = np.dot(coef, sigma) / np.dot(coef, coef)
    return energy