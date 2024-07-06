import time
import numpy as np
from fcipy.lib import ci
from fcipy.printing import print_dav_iteration_header, print_dav_iteration, dav_calculation_summary, get_timestamp

def davidson_ci(system, det, e1int, e2int, b, e, e_diag, maxit=200, max_size=30, tolerance=1.0e-08):

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
            print_dav_iteration(niter, e + system.nuclear_repulsion, resnorm, delta_e, elapsed_time)
            is_converged = True
            break
        # DPR residual update
        for i in range(ndim):
            r[i] /= (e_diag[i] - e + 1.0e-09)
        r /= np.linalg.norm(r)
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
        print_dav_iteration(niter, e + system.nuclear_repulsion, resnorm, delta_e, elapsed_time)
        curr_size += 1

    return e + system.nuclear_repulsion, v, is_converged

# def davidson_dense(A, b, e, maxit=200, max_size=30, tolerance=1.0e-08):
#
#     B = np.zeros((len(b), max_size))
#     S = np.zeros((len(b), max_size))
#     G = np.zeros((max_size, max_size))
#
#     B[:, 0] = b
#     S[:, 0] = np.dot(A, b)
#
#     is_converged = False
#     curr_size = 1
#     for niter in range(maxit):
#         # start iteration
#         t1 = time.time()
#         # store old energy
#         e_old = e.copy()
#         # Form and diagonalize subspace matrix
#         for p in range(curr_size):
#             G[curr_size - 1, p] = np.dot(B[:, curr_size - 1].T, S[:, p])
#             G[p, curr_size - 1] = np.dot(S[:, curr_size - 1].T, B[:, p])
#         e, alpha = np.linalg.eig(G[:curr_size, :curr_size])
#         # Pick the eigenvalue
#         idx = np.argsort(abs(alpha[0, :]))
#         alpha = np.real(alpha[:, idx[-1]])
#         e = np.real(e[idx[-1]])
#         # Compute eigenvector and residual
#         v = np.dot(B[:, :curr_size], alpha)
#         r = np.dot(S[:, :curr_size], alpha) - e * v
#         # residual and convergence measures
#         resnorm = np.linalg.norm(r)
#         delta_e = e - e_old
#         # check for convergence and exit
#         if resnorm < tolerance and abs(delta_e) < tolerance:
#             is_converged = True
#             elapsed_time = time.time() - t1
#             print("iteration: %d   residual: %.9f   eigenvalue:  %.9f   elapsed time:  %.2f" % (niter, resnorm, e, elapsed_time))
#             print("converged!")
#             break
#         # DPR residual update
#         for i in range(A.shape[0]):
#             r[i] /= (A[i, i] - e + 1.0e-09)
#         r /= np.linalg.norm(r)
#         for p in range(curr_size):
#             b = B[:, p] / np.linalg.norm(B[:, p])
#             r -= np.dot(b.T, r) * b
#         r /= np.linalg.norm(r)
#         # enlarge the subsapce
#         if curr_size < max_size:
#             B[:, curr_size] = r
#             S[:, curr_size] = np.dot(A, r)
#         else: # if max_size exceeded, restart from current guess to eigenvector
#             B[:, 0] = v
#             S[:, 0] = np.dot(A, v)
#             curr_size = 0
#         # update iteration counter
#         elapsed_time = time.time() - t1
#         print("iteration: %d   residual: %.9f   eigenvalue:  %.9f   elapsed time:  %.2f s" % (niter, resnorm, e, elapsed_time))
#         curr_size += 1
#     else:
#         print("root not converged")
#
#     return e, v
#
# def run_davidson_dense(A, nroot):
#
#     ndim = A.shape[0]
#     A_diag = np.diag(A)
#     idx = np.argsort(A_diag)
#
#     e = np.zeros(nroot)
#     v = np.zeros((ndim, nroot))
#
#     for n in range(nroot):
#         print("Solving for root", n + 1)
#         e0 = A_diag[idx[n]]
#         b0 = np.zeros(ndim)
#         b0[idx[n]] = 1.0
#
#         print("start eigenvalue:", e0)
#         e[n], v[:, n] = davidson(A, b0, e0)
#     print("completed all Davidson iterations")
#
#     return e, v

def run_davidson(system, det, e1int, e2int, nroot):
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
        print("\n   Energy of initial guess = {:>10.10f}".format(e0 + system.nuclear_repulsion))
        e[n], v[:, n], is_converged = davidson_ci(system, det, e1int, e2int, b0, e0, A_diag)
        dav_calculation_summary(e[n], is_converged, n, system, 0.09)

    return e, v

