"""Main calculation driver module of FCIpy."""

import numpy as np
import scipy
from fcipy.interfaces import load_pyscf_integrals, load_gamess_integrals, general_system

class Driver:

    @classmethod
    def from_custom(cls, nelectrons, norbitals, nfrozen, ndelete, mult, e1int, e2int, nuclear_repulsion, point_group, orbital_symmetry):
        return cls(
                    *general_system(nelectrons, norbitals, nfrozen, ndelete, mult, e1int, e2int, nuclear_repulsion, point_group, orbital_symmetry)
                  )

    @classmethod
    def from_pyscf(cls, meanfield, nfrozen, ndelete=0):
        return cls(
                    *load_pyscf_integrals(meanfield, nfrozen, ndelete)
                  )

    @classmethod
    def from_gamess(cls, logfile, nfrozen, ndelete=0, multiplicity=None, fcidump=None, onebody=None, twobody=None):
        return cls(
                    *load_gamess_integrals(logfile, fcidump, onebody, twobody, nfrozen, ndelete, multiplicity)
                   )

    def __init__(self, system, e1int, e2int):
        self.system = system
        self.e1int = e1int
        self.e2int = e2int
        self.coef = None
        self.coef_left = None
        self.total_energy = None
        self.state_eigval = None
        self.det = None
        self.ndet = 0
        self.Hmat = None
        self.N_int = int(np.floor(self.system.norbitals / 64) + 1)
        self.num_alpha = 0
        self.num_beta = 0
        self.rdm1 = [None for _ in range(100)]
        self.rdms = [{} for _ in range(100)]
        self.nat_orb = [None for _ in range(100)]
        self.nat_occ_num = [None for _ in range(100)]

    def load_determinants(self, file=None, max_excit_rank=-1, target_irrep=None):
        from fcipy.determinants import read_determinants
        from fcipy.determinants import fci_space
        if file is not None:
            self.det = read_determinants(file)
        else:
            self.det = fci_space(self.system, max_excit_rank=max_excit_rank, target_irrep=target_irrep)

        assert self.N_int == self.det.shape[0]
        self.ndet = self.det.shape[2]

        # record number of unique alpha and beta strings
        self.num_alpha = len(np.unique(self.det[:, 0, :]))
        self.num_beta = len(np.unique(self.det[:, 1, :]))

    def print_determinants(self):
        for idet in range(self.ndet):
            print(f"determinant {idet + 1}: alpha = {self.det[:, 0, idet]}  beta = {self.det[:, 1, idet]}")

    def print_ci_vector(self, state=0, prtol=0.01, file=None):
        from fcipy.printing import print_ci_amplitudes, print_ci_amplitudes_to_file
        if file is None:
            print("   Right CI state")
            print_ci_amplitudes(self.system, self.det, self.coef[:, state], thresh=prtol)
            if self.coef_left[:, state] is not None:
                print("   Left CI state")
                print_ci_amplitudes(self.system, self.det, self.coef_left[:, state], thresh=prtol)
        else:
            print_ci_amplitudes_to_file(file, self.system, self.det, self.coef[:, state], thresh=prtol)

    def run_ci(self, nroot, convergence=1.0e-08, max_size=30, maxit=200, opt=True, herm=True, prtol=0.09):
        from fcipy.davidson import run_davidson, run_davidson_opt
        if opt:
            # for now, the sorting routine in the optimized CI requires that N_int = 1
            assert self.N_int == 1
            self.total_energy, self.state_eigval, self.coef = run_davidson_opt(self.system, self.det, self.num_alpha, self.num_beta, self.e1int, self.e2int, self.coef,
                                                            nroot, convergence=convergence, max_size=max_size, maxit=maxit, herm=herm, print_thresh=prtol)
        else:
            self.total_energy, self.state_eigval, self.coef = run_davidson(self.system, self.det, self.e1int, self.e2int, self.coef, nroot,
                                                        convergence=convergence, max_size=max_size, maxit=maxit, herm=herm, print_thresh=prtol)

    def build_hamiltonian(self, opt=True, herm=True):
        from fcipy.hamiltonian import build_hamiltonian
        self.Hmat = build_hamiltonian(self.det, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta, herm=herm)

    def diagonalize_hamiltonian(self, herm=True):

        if herm:
            self.state_eigval, self.coef = scipy.linalg.eigh(self.Hmat)
            self.coef_left = self.coef.copy()
        else:
            self.state_eigval, self.coef_left, self.coef = scipy.linalg.eig(self.Hmat, right=True, left=True)
            idx = np.argsort(self.state_eigval)
            self.state_eigval = self.state_eigval[idx]
            self.coef = self.coef[:, idx]
            self.coef_left = self.coef_left[:, idx]

            print(f"   [fcipy]: norm of imag(C_R) = {np.linalg.norm(np.imag(self.coef))}")
            print(f"   [fcipy]: norm of imag(C_L) = {np.linalg.norm(np.imag(self.coef_left))}")

            # Biorthonormalize L and R manually
            # M = np.conj(self.coef_left).T @ self.coef
            # M_L, M_U = scipy.linalg.lu(M, permute_l=True, overwrite_a=True, check_finite=False)
            # self.coef_left = scipy.linalg.inv(M_L, overwrite_a=True, check_finite=False) @ self.coef_left
            # self.coef = self.coef @ scipy.linalg.inv(M_U, overwrite_a=True, check_finite=False)

            # Check biorthonormality
            # LR = np.dot(np.conj(self.coef_left).T, self.coef)
            # assert np.allclose(LR, np.eye(self.coef.shape[1]), atol=1.0e-09, rtol=1.0e-09), "Left- and right-eigenvectors of non-Hermitian Hamiltonian are not biorthonormal!"

        self.total_energy = self.state_eigval + self.system.frozen_energy + self.system.nuclear_repulsion

    def compute_energy(self, herm):
        from fcipy.ci_energy import energy_ci_state

        energy = np.zeros(self.coef.shape[1])
        for istate in range(self.coef.shape[1]):
            energy[istate] = energy_ci_state(self.det, self.coef[:, istate], self.e1int, self.e2int,
                                             self.system.noccupied_alpha, self.system.noccupied_beta,
                                             self.num_alpha, self.num_beta, herm)
        # Add frozen energy and nuclear repulsion to eigenvalues to make total energy
        total_energy = energy + self.system.frozen_energy + self.system.nuclear_repulsion
        return energy, total_energy

    def apply_hamiltonian(self, c, herm):
        from fcipy.lib import ci

        e_diag = ci.ci.calc_diagonal(self.det, self.system.noccupied_alpha, self.system.noccupied_beta, self.e1int, self.e2int)
        if herm:
            sigma = ci.ci.calc_sigma_opt(self.det, e_diag, c, self.num_alpha, self.num_beta, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta)
        else:
            sigma = ci.ci.calc_sigma_nonhermitian_opt(self.det, e_diag, c, self.num_alpha, self.num_beta, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta)

        return sigma

    def build_rdm1s(self, i=0):
        from fcipy.density import compute_rdm1s

        if self.coef_left is None:
            vL = self.coef[:, i]
        else:
            vL = self.coef_left[:, i].conj()

        dm1a, dm1b = compute_rdm1s(self.det, self.coef[:, i], vL, self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
        self.rdms[i]['a'] = dm1a
        self.rdms[i]['b'] = dm1b

    def build_rdm2s(self, i=0):
        from fcipy.density import compute_rdm2s

        if self.coef_left is None:
            vL = self.coef[:, i]
        else:
            vL = self.coef_left[:, i].conj()

        dm2aa, dm2ab, dm2bb = compute_rdm2s(self.det, self.coef[:, i], vL, self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
        self.rdms[i]['aa'] = dm2aa
        self.rdms[i]['ab'] = dm2ab
        self.rdms[i]['bb'] = dm2bb

    def build_rdm3s(self, i):
        from fcipy.density import compute_rdm3s

        if self.coef_left is None:
            vL = self.coef[:, i]
        else:
            vL = self.coef_left[:, i].conj()

        dm3aaa, dm3aab, dm3abb, dm3bbb = compute_rdm3s(self.det, self.coef[:, i], vL, self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
        self.rdms[i]['aaa'] = dm3aaa
        self.rdms[i]['aab'] = dm3aab
        self.rdms[i]['abb'] = dm3abb
        self.rdms[i]['bbb'] = dm3bbb

    def compute_rdm123s(self, i, use_right=False):
        from fcipy.density import compute_rdm1s, compute_rdm2s, compute_rdm3s
        rdms_i = {}

        if use_right:
            vL = self.coef[:, i].copy()
        else:
            if self.coef_left is None:
                vL = self.coef[:, i]
            else:
                vL = self.coef_left[:, i].conj()

        rdms_i['a'], rdms_i['b'] = compute_rdm1s(self.det, self.coef[:, i], vL, self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
        rdms_i['aa'], rdms_i['ab'], rdms_i['bb'] = compute_rdm2s(self.det, self.coef[:, i], vL, self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
        rdms_i['aaa'], rdms_i['aab'], rdms_i['abb'], rdms_i['bbb'] = compute_rdm3s(self.det, self.coef[:, i], vL, self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)

        return rdms_i

