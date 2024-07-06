import numpy as np

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

def fci(system, target_irrep):
    from itertools import combinations

    sym_target = system.point_group_irrep_to_number[target_irrep]

    N_int = int(np.floor(system.norbitals / 64) + 1)
    ndet = get_fci_space_dimension(system, target_irrep)

    print("   Dimension of FCI space = ", ndet)

    det = np.zeros((N_int, 2, ndet), dtype=np.int64)

    orbs = np.arange(system.norbitals)

    kout = 0
    for alpha in combinations(orbs, system.noccupied_alpha):

        # sym_a = 0
        # for i, a in enumerate(alpha):
        #     sym = system.point_group_irrep_to_number[system.orbital_symmetries[a]]
        #     sym_a = sym_a ^ sym
        # if sym_a != sym_target: continue

        I_a = sum([2**n for n in alpha])

        for beta in combinations(orbs, system.noccupied_beta):

            # sym_b = 0
            # for i, b in enumerate(beta):
            #     sym = system.point_group_irrep_to_number[system.orbital_symmetries[b]]
            #     sym_b = sym_b ^ sym
            # if sym_b != sym_target: continue

            I_b = sum([2**n for n in beta])

            det[0, 0, kout] = I_a
            det[0, 1, kout] = I_b
            kout += 1

    return det

def get_fci_space_dimension(system, target_irrep):
    from itertools import combinations

    orbs = np.arange(system.norbitals)
    sym_target = system.point_group_irrep_to_number[target_irrep]

    kout = 0
    for alpha in combinations(orbs, system.noccupied_alpha):

        # sym_a = 0
        # for i, a in enumerate(alpha):
        #     sym = system.point_group_irrep_to_number[system.orbital_symmetries[a]]
        #     sym_a = sym_a ^ sym
        # if sym_a != sym_target: continue

        for beta in combinations(orbs, system.noccupied_beta):

            # sym_b = 0
            # for i, b in enumerate(beta):
            #     sym = system.point_group_irrep_to_number[system.orbital_symmetries[b]]
            #     sym_b = sym_b ^ sym
            # if sym_b != sym_target: continue

            kout += 1

    return kout
