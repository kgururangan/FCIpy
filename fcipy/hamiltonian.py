import numpy as np
from fcipy.lib import ci

def build_hamiltonian(det, e1int, e2int, noa, nob):
    Hmat = ci.ci.build_hamiltonian(det, e1int, e2int, noa, nob)
    return Hmat