import time
import numpy as np
from fcipy.lib import ci
from fcipy.printing import print_dav_iteration_header, print_dav_iteration, dav_calculation_summary, get_timestamp

def davidson(system, det, e1int, e2int, b, e, e_diag, tolerance, max_size, maxit):

    print_dav_iteration_header()

    ndim = det.shape[2]
    noa = system.noccupied_alpha
    nob = system.noccupied_beta

    B = np.zeros((len(b), max_size))
    S = np.zeros((len(b), max_size))
    G = np.zeros((max_size, max_size))

    B[:, 0] = b
    S[:, 0] = ci.ci.calc_sigma(det, e1int, e2int, noa, nob, b)

    is_converged = False
    curr_size = 1
    for niter in range(maxit):
        # start iteration
        t1 = time.time()
        # store old energy
        e_old = e.copy()
        # Form and diagonalize subspace matrix
        for p in range(curr_size):
            G[curr_size - 1, p] = np.dot(B[:, curr_size - 1].T, S[:, p])
            G[p, curr_size - 1] = np.dot(S[:, curr_size - 1].T, B[:, p])
        e, alpha = np.linalg.eig(G[:curr_size, :curr_size])
        # Pick the eigenvalue
        idx = np.argsort(abs(alpha[0, :]))
        alpha = np.real(alpha[:, idx[-1]])
        e = np.real(e[idx[-1]])
        # Compute eigenvector and residual
        v = np.dot(B[:, :curr_size], alpha)
        r = np.dot(S[:, :curr_size], alpha) - e * v
        # residual and convergence measures
        resnorm = np.linalg.norm(r)
        delta_e = e - e_old
        # check for convergence and exit
        if resnorm < tolerance and abs(delta_e) < tolerance:
            is_converged = True
            elapsed_time = time.time() - t1
            print_dav_iteration(niter, e + system.nuclear_repulsion + system.frozen_energy, resnorm, delta_e, elapsed_time)
            is_converged = True
            break
        # DPR residual update
        r /= (e_diag - e + 1.0e-012)
        # Orthogonalize residual against trial vectors
        for p in range(curr_size):
            b = B[:, p] / np.linalg.norm(B[:, p])
            r -= np.dot(b.T, r) * b
        r /= np.linalg.norm(r)
        # enlarge the subsapce
        if curr_size < max_size:
            B[:, curr_size] = r
            S[:, curr_size] = ci.ci.calc_sigma(det, e1int, e2int, noa, nob, r) 
        else: # if max_size exceeded, restart from current guess to eigenvector
            B[:, 0] = v
            S[:, 0] = ci.ci.calc_sigma(det, e1int, e2int, noa, nob, v) 
            curr_size = 0
        # update iteration counter
        elapsed_time = time.time() - t1
        print_dav_iteration(niter, e + system.nuclear_repulsion + system.frozen_energy, resnorm, delta_e, elapsed_time)
        curr_size += 1

    return e + system.nuclear_repulsion + system.frozen_energy, v, is_converged

def davidson_opt(system, det, num_alpha, num_beta, e1int, e2int, b, e, e_diag, tolerance, max_size, maxit):

    print_dav_iteration_header()

    ndim = det.shape[2]
    noa = system.noccupied_alpha
    nob = system.noccupied_beta

    B = np.zeros((len(b), max_size))
    S = np.zeros((len(b), max_size))
    G = np.zeros((max_size, max_size))

    B[:, 0] = b
    det, B[:, 0], S[:, 0], e_diag = ci.ci.calc_sigma_opt(det, num_alpha, num_beta, e1int, e2int, noa, nob, b, e_diag)

    is_converged = False
    curr_size = 1
    for niter in range(maxit):
        # start iteration
        t1 = time.time()
        # store old energy
        e_old = e.copy()
        # Form and diagonalize subspace matrix
        for p in range(curr_size):
            G[curr_size - 1, p] = np.dot(B[:, curr_size - 1].T, S[:, p])
            G[p, curr_size - 1] = np.dot(S[:, curr_size - 1].T, B[:, p])
        e, alpha = np.linalg.eig(G[:curr_size, :curr_size])
        # Pick the eigenvalue
        idx = np.argsort(abs(alpha[0, :]))
        alpha = np.real(alpha[:, idx[-1]])
        e = np.real(e[idx[-1]])
        # Compute eigenvector and residual
        v = np.dot(B[:, :curr_size], alpha)
        r = np.dot(S[:, :curr_size], alpha) - e * v
        # residual and convergence measures
        resnorm = np.linalg.norm(r)
        delta_e = e - e_old
        # check for convergence and exit
        if resnorm < tolerance and abs(delta_e) < tolerance:
            elapsed_time = time.time() - t1
            print_dav_iteration(niter, e + system.nuclear_repulsion + system.frozen_energy, resnorm, delta_e, elapsed_time)
            is_converged = True
            break
        # DPR residual update
        r /= (e_diag - e + 1.0e-012)
        # Orthogonalize residual against trial vectors
        for p in range(curr_size):
            b = B[:, p] / np.linalg.norm(B[:, p])
            r -= np.dot(b.T, r) * b
        r /= np.linalg.norm(r)
        # enlarge the subsapce
        if curr_size < max_size:
            B[:, curr_size] = r
            det, B[:, curr_size], S[:, curr_size], e_diag = ci.ci.calc_sigma_opt(det, num_alpha, num_beta, e1int, e2int, noa, nob, r, e_diag)
        else: # if max_size exceeded, restart from current guess to eigenvector
            B[:, 0] = v
            det, B[:, 0], S[:, 0], e_diag = ci.ci.calc_sigma_opt(det, num_alpha, num_beta, e1int, e2int, noa, nob, v, e_diag)
            curr_size = 0
        # update iteration counter
        elapsed_time = time.time() - t1
        print_dav_iteration(niter, e + system.nuclear_repulsion + system.frozen_energy, resnorm, delta_e, elapsed_time)
        curr_size += 1

    return det, e + system.nuclear_repulsion + system.frozen_energy, v, is_converged

# def block_davidson(system, det, e1int, e2int, b, e, e_diag, tolerance, max_size, maxit):
#
#     print_dav_iteration_header()
#
#     # Number of roots
#     nroot = b.shape[1]
#     ndim = b.shape[0]
#
#     # Allocate the B (correction/subspace), sigma (HR), and G (interaction) matrices
#     sigma = np.zeros((ndim, max_size))
#     B = np.zeros((ndim, max_size))
#     G = np.zeros((max_size, max_size))
#
#     # Initial values
#     num_add = 0
#     curr_size = 0
#     for j, istate in enumerate(state_index):
#         B[:, j] = b[:, j]
#         S[:, j] = ci.ci.calc_sigma(det, e1int, e2int, noa, nob, b[:, j])
#         num_add += 1
#         curr_size += 1
#
#     is_converged = [False] * nroot
#     residual = np.zeros(nroot)
#     delta_energy = np.zeros(nroot)
#     for niter in range(maxit):
#         # start iteration
#         t1 = time.time()
#         # store old energy
#         e_old = e.copy()
#         # Form and diagonalize subspace matrix
#         for p in range(curr_size):
#             for j in range(1, num_add + 1):
#                 G[curr_size - j, p] = np.dot(B[:, curr_size - j].T, S[:, p])
#                 G[p, curr_size - j] = np.dot(S[:, curr_size - j].T, B[:, p])
#         e, alpha_full = np.linalg.eig(G[:curr_size, :curr_size])
#         #
#         num_add = 0
#         nmax_add = sum([not x for x in is_converged])
#         alpha = np.zeros((curr_size, nroot))
#         for j in range(nroot):
#             if is_converged[j]: continue
#             idx = np.argsort(abs(alpha_full[j, :]))
#             iselect = idx[-1]
#             alpha[:, j] = np.real(alpha_full[:, iselect])
#             e[j] = np.real(e[iselect])
#             # Compute eigenvector and residual
#             v[:, j] = np.dot(B[:, :curr_size], alpha)
#             r = np.dot(S[:, :curr_size], alpha) - e * v
#             # residual and convergence measures
#             resnorm = np.linalg.norm(r)
#             delta_e = e - e_old
#         # check for convergence and exit
#         if resnorm < tolerance and abs(delta_e) < tolerance:
#             is_converged = True
#             elapsed_time = time.time() - t1
#             print_dav_iteration(niter, e + system.nuclear_repulsion, resnorm, delta_e, elapsed_time)
#             is_converged = True
#             break
#         # DPR residual update
#         for i in range(ndim):
#             r[i] /= (e_diag[i] - e + 1.0e-09)
#         r /= np.linalg.norm(r)
#         for p in range(curr_size):
#             b = B[:, p] / np.linalg.norm(B[:, p])
#             r -= np.dot(b.T, r) * b
#         r /= np.linalg.norm(r)
#         # enlarge the subsapce
#         if curr_size < max_size:
#             B[:, curr_size] = r
#             S[:, curr_size] = ci.ci.calc_sigma(det, e1int, e2int, noa, nob, r)
#         else: # if max_size exceeded, restart from current guess to eigenvector
#             B[:, 0] = v
#             S[:, 0] = ci.ci.calc_sigma(det, e1int, e2int, noa, nob, v)
#             curr_size = 0
#         # update iteration counter
#         elapsed_time = time.time() - t1
#         print_dav_iteration(niter, e + system.nuclear_repulsion, resnorm, delta_e, elapsed_time)
#         curr_size += 1
#     return e + system.nuclear_repulsion, v, is_converged
#
#             r = np.dot(B[:, :curr_size], alpha[:, j])
#
#             # calculate residual vector: r_i = S_{iK}*alpha_{K} - omega * r_i
#             R[istate].unflatten(np.dot(sigma[:, :curr_size], alpha[:, j]) - omega[istate] * r)
#             residual[j] = np.linalg.norm(R[istate].flatten())
#             delta_energy[j] = omega[istate] - omega_old[istate]
#
#             # Check convergence
#             if residual[j] < options["amp_convergence"] and abs(delta_energy[j]) < options["energy_convergence"]:
#                 is_converged[j] = True
#             else:
#                 # update the residual vector
#                 if t3_excitations and r3_excitations:
#                     R[istate] = update_r(R[istate], omega[istate], H, options["RHF_symmetry"], system, r3_excitations)
#                 else:
#                     R[istate] = update_r(R[istate], omega[istate], H, options["RHF_symmetry"], system)
#                 q = R[istate].flatten()
#                 for p in range(curr_size + num_add):
#                     b = B[:, p] / np.linalg.norm(B[:, p])
#                     q -= np.dot(b.T, q) * b
#                 q /= np.linalg.norm(q)
#                 R[istate].unflatten(q)
#                 B[:, curr_size + num_add] = q
#                 sigma[:, curr_size + num_add] = HR(dR, R[istate], T, H, options["RHF_symmetry"], system)
#                 num_add += 1
#
#             # Store the root you've solved for
#             R[istate].unflatten(r)
#
#         # Check for all roots converged and break
#         if all(is_converged):
#             print("   All roots converged")
#             break
#
#         # print the iteration
#         elapsed_time = time.perf_counter() - t1
#         print_block_eomcc_iteration(niter + 1, curr_size, omega, residual, delta_energy, elapsed_time, state_index)
#         curr_size += num_add
#
#     return R, omega, is_converged

def run_davidson(system, det, e1int, e2int, nroot, convergence, max_size, maxit, print_thresh=0.09):

    ndim = det.shape[2]
    A_diag = ci.ci.calc_diagonal(det, system.noccupied_alpha, system.noccupied_beta, e1int, e2int)
    idx = np.argsort(A_diag)

    e = np.zeros(nroot)
    v = np.zeros((ndim, nroot))

    for n in range(nroot):
        e0 = A_diag[idx[n]]
        b0 = np.zeros(ndim)
        b0[idx[n]] = 1.0

        print("   CI calculation for root %d started on" % n, get_timestamp())
        print("\n   Energy of initial guess = {:>10.10f}".format(e0 + system.nuclear_repulsion + system.frozen_energy))
        e[n], v[:, n], is_converged = davidson(system, det, e1int, e2int, b0, e0, A_diag, convergence, max_size, maxit)
        dav_calculation_summary(e[n], det, v[:, n], is_converged, n, system, print_thresh)
    return e, v

def run_davidson_opt(system, det, num_alpha, num_beta, e1int, e2int, nroot, convergence, max_size, maxit, print_thresh=0.09):

    ndim = det.shape[2]
    e = np.zeros(nroot)
    v = np.zeros((ndim, nroot))

    for n in range(nroot):
        A_diag = ci.ci.calc_diagonal(det, system.noccupied_alpha, system.noccupied_beta, e1int, e2int)
        idx = np.argsort(A_diag)
        e0 = A_diag[idx[n]]
        b0 = np.zeros(ndim)
        b0[idx[n]] = 1.0

        print("   CI calculation for root %d started on" % n, get_timestamp())
        print(f"   Number of unique alpha strings: {num_alpha}")
        print(f"   Number of unique beta strings: {num_beta}")
        print("\n   Energy of initial guess = {:>10.10f}".format(e0 + system.nuclear_repulsion + system.frozen_energy))
        det, e[n], v[:, n], is_converged = davidson_opt(system, det, num_alpha, num_beta, e1int, e2int, b0, e0, A_diag, convergence, max_size, maxit)
        dav_calculation_summary(e[n], det, v[:, n], is_converged, n, system, print_thresh)
    return det, e, v
