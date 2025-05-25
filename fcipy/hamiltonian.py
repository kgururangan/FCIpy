from fcipy.lib import ci

def build_hamiltonian(det, e1int, e2int, noa, nob, herm):
    if herm:
        Hmat = ci.ci.build_hamiltonian(det, e1int, e2int, noa, nob)
    else:
        Hmat = ci.ci.build_hamiltonian_nonhermitian(det, e1int, e2int, noa, nob)
    return Hmat
