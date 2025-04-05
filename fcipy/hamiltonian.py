from fcipy.lib import ci

def build_hamiltonian(det, e1int, e2int, noa, nob):
    Hmat = ci.ci.build_hamiltonian(det, e1int, e2int, noa, nob)
    return Hmat

# def build_hamiltonian_opt(det, num_alpha, num_beta, e1int, e2int, noa, nob, e_ref):
#     det, Hmat = ci.ci.build_hamiltonian_opt(det, num_alpha, num_beta, e1int, e2int, noa, nob, e_ref)
#     return det, Hmat
