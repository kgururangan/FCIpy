import numpy as np
from fcipy.lib import ci

def compute_rdm1s(det, coef, coef_left, mo_num, noa, nob):
    return ci.ci.compute_rdm1s(det, coef, coef_left, mo_num, noa, nob)

def compute_rdm2s(det, coef, coef_left, mo_num, noa, nob):
    return ci.ci.compute_rdm2s(det, coef, coef_left, mo_num, noa, nob)

def compute_rdm3s(det, coef, coef_left, mo_num, noa, nob):
    return ci.ci.compute_rdm3s(det, coef, coef_left, mo_num, noa, nob)
