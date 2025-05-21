import numpy as np
from fcipy.lib import ci

def compute_rdm1(det, coef, mo_num):
    return ci.ci.compute_rdm1(det, coef, mo_num)

def compute_rdm1s(det, coef, mo_num, noa, nob):
    return ci.ci.compute_rdm1s(det, coef, mo_num, noa, nob)

def compute_rdm2s(det, coef, mo_num, noa, nob):
    return ci.ci.compute_rdm2s(det, coef, mo_num, noa, nob)
