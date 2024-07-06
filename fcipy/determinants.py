import time
import numpy as np
from fcipy.excitations import get_excitation_degree

def read_determinants(file):
    with open(file, "r") as f:
        for i, line in enumerate(f.readlines()):
            if i == 0: # first line stores N_int, something, ndet
                N_int, _, ndet = [int(x) for x in line.split()]
                det = np.zeros(N_int * 2 * ndet, dtype=np.int64, order="F")
                continue
            # load the numbers into a flat array
            det[i - 1] = np.int64(line.split()[0])

    # Reshape the array in Fortran order (F order is important!)
    det = np.reshape(det, (N_int, 2, ndet), order="F")
    return det

def check_symmetry(sym_alpha, sym_beta, sym_target):
    if sym_target == -1:
        return True
    else:
        return sym_alpha ^ sym_beta == sym_target

def fci_space(system, max_excit_rank, target_irrep):
    from itertools import combinations

    t_start = time.perf_counter()

    if target_irrep is None:
        sym_target = -1
    else:
        sym_target = system.point_group_irrep_to_number[target_irrep]
    if max_excit_rank == -1:
        max_excit_rank = system.nelectrons

    print("   Generating the CI space")
    print("   ------------------------")
    print("   Number of occupied alpha = ", system.noccupied_alpha)
    print("   Number of occupied beta = ", system.noccupied_beta)
    print(f"   Target Irrep = {target_irrep} ({system.point_group})")
    print(f"   Maximum Excitation Rank = {max_excit_rank}")

    N_int = int(np.floor(system.norbitals / 64) + 1)
    ndet = get_fci_space_dimension(system, max_excit_rank, target_irrep)

    # Hartree-Fock determinant
    I_a_ref = np.zeros(N_int, dtype=np.int64)
    for n in range(system.noccupied_alpha):
        I_a_ref[n // 64] += 2 ** (n % 64)
    I_b_ref = np.zeros(N_int, dtype=np.int64)
    for n in range(system.noccupied_beta):
        I_b_ref[n // 64] += 2 ** (n % 64)

    print("   Reference determinant = ", (I_a_ref, I_b_ref))
    print("   Dimension of CI space = ", ndet)

    det = np.zeros((N_int, 2, ndet), dtype=np.int64)
    orbs = np.arange(system.norbitals)
    kout = 0
    for alpha in combinations(orbs, system.noccupied_alpha):

        sym_a = 0
        for i, a in enumerate(alpha):
            sym = system.point_group_irrep_to_number[system.orbital_symmetries[a]]
            sym_a = sym_a ^ sym

        I_a = np.zeros(N_int, dtype=np.int64)
        for n in alpha:
            I_a[n // 64] += 2**(n % 64)

        degree_a = get_excitation_degree(I_a, I_a_ref)

        for beta in combinations(orbs, system.noccupied_beta):

            sym_b = 0
            for i, b in enumerate(beta):
                sym = system.point_group_irrep_to_number[system.orbital_symmetries[b]]
                sym_b = sym_b ^ sym

            I_b = np.zeros(N_int, dtype=np.int64)
            for n in beta:
                I_b[n // 64] += 2 ** (n % 64)
            degree_b = get_excitation_degree(I_b, I_b_ref)

            degree = degree_a + degree_b

            if check_symmetry(sym_a, sym_b, sym_target) and degree <= max_excit_rank:
                det[:, 0, kout] = I_a
                det[:, 1, kout] = I_b
                kout += 1
    t_end = time.perf_counter()
    print(f"   Completed in {t_end - t_start} seconds\n")
    return det

def get_fci_space_dimension(system, max_excit_rank, target_irrep):
    from itertools import combinations

    if target_irrep is None:
        sym_target = -1
    else:
        sym_target = system.point_group_irrep_to_number[target_irrep]

    if max_excit_rank == -1:
        max_excit_rank = system.nelectrons

    N_int = int(np.floor(system.norbitals / 64) + 1)

    # Hartree-Fock determinant
    I_a_ref = np.zeros(N_int, dtype=np.int64)
    for n in range(system.noccupied_alpha):
        I_a_ref[n // 64] += 2 ** (n % 64)
    I_b_ref = np.zeros(N_int, dtype=np.int64)
    for n in range(system.noccupied_beta):
        I_b_ref[n // 64] += 2 ** (n % 64)

    orbs = np.arange(system.norbitals)
    kout = 0
    for alpha in combinations(orbs, system.noccupied_alpha):

        sym_a = 0
        for i, a in enumerate(alpha):
            sym = system.point_group_irrep_to_number[system.orbital_symmetries[a]]
            sym_a = sym_a ^ sym

        I_a = np.zeros(N_int, dtype=np.int64)
        for n in alpha:
            I_a[n // 64] += 2**(n % 64)

        degree_a = get_excitation_degree(I_a, I_a_ref)

        for beta in combinations(orbs, system.noccupied_beta):

            sym_b = 0
            for i, b in enumerate(beta):
                sym = system.point_group_irrep_to_number[system.orbital_symmetries[b]]
                sym_b = sym_b ^ sym

            I_b = np.zeros(N_int, dtype=np.int64)
            for n in beta:
                I_b[n // 64] += 2 ** (n % 64)
            degree_b = get_excitation_degree(I_b, I_b_ref)

            degree = degree_a + degree_b

            if check_symmetry(sym_a, sym_b, sym_target) and degree <= max_excit_rank:
                kout += 1
    return kout
