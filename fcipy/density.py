import numpy as np
from fcipy.lib import ci

def compute_rdm1(det, coef, mo_num):
    return ci.ci.compute_density_matrix(det, coef, mo_num)
