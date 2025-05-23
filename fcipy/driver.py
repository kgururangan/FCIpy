"""Main calculation driver module of FCIpy."""

import numpy as np
from fcipy.interfaces import load_pyscf_integrals, load_gamess_integrals

class Driver:

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

    def __init__(self, system, e1int, e2int, herm=True):
        self.system = system
        self.herm = herm
        self.e1int = e1int
        self.e2int = e2int
        self.coef = None
        self.total_energy = None
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
            print_ci_amplitudes(self.system, self.det, self.coef[:, state], thresh=prtol)
        else:
            print_ci_amplitudes_to_file(file, self.system, self.det, self.coef[:, state], thresh=prtol)

    def run_ci(self, nroot, convergence=1.0e-08, max_size=30, maxit=200, opt=True, prtol=0.09):
        from fcipy.davidson import run_davidson, run_davidson_opt
        if opt:
            # for now, the sorting routine in the optimized CI requires that N_int = 1
            assert self.N_int == 1
            self.det, self.total_energy, self.coef = run_davidson_opt(self.system, self.det, self.num_alpha, self.num_beta, self.e1int, self.e2int, nroot,
                                                                      convergence=convergence, max_size=max_size, maxit=maxit, print_thresh=prtol, herm=self.herm)
        else:
            self.total_energy, self.coef = run_davidson(self.system, self.det, self.e1int, self.e2int, nroot,
                                                        convergence=convergence, max_size=max_size, maxit=maxit, print_thresh=prtol, herm=self.herm)

    def build_hamiltonian(self, opt=True):
        from fcipy.hamiltonian import build_hamiltonian
        # if opt:
        #     # for now, the sorting routine in the optimized CI requires that N_int = 1
        #     assert self.N_int == 1
        #     print(self.det.shape)
        #     self.det, self.Hmat = build_hamiltonian_opt(self.det, self.num_alpha, self.num_beta, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta, self.system.reference_energy)
        # else:
        self.Hmat = build_hamiltonian(self.det, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta, herm=self.herm)

    def diagonalize_hamiltonian(self, opt=True):

        if self.Hmat is None:
            self.build_hamiltonian(opt=opt)

        if self.herm:
            self.total_energy, self.coef = np.linalg.eigh(self.Hmat)
        else:
            self.total_energy, self.coef = np.linalg.eig(self.Hmat)
            idx = np.argsort(self.total_energy)
            self.total_energy = self.total_energy[idx]
            self.coef = self.coef[:, idx]

        self.total_energy += self.system.frozen_energy
        self.total_energy += self.system.nuclear_repulsion

    def build_rdm1s(self):
        from fcipy.density import compute_rdm1s
        for i in range(self.coef.shape[1]):
            dm1a, dm1b = compute_rdm1s(self.det, self.coef[:, i], self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
            self.rdms[i]['a'] = dm1a
            self.rdms[i]['b'] = dm1b

    def build_rdm2s(self):
        from fcipy.density import compute_rdm2s
        for i in range(self.coef.shape[1]):
            dm2aa, dm2ab, dm2bb = compute_rdm2s(self.det, self.coef[:, i], self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
            self.rdms[i]['aa'] = dm2aa
            self.rdms[i]['ab'] = dm2ab
            self.rdms[i]['bb'] = dm2bb

    def build_rdm3s(self):
        from fcipy.density import compute_rdm3s
        for i in range(self.coef.shape[1]):
            dm3aaa, dm3aab, dm3abb, dm3bbb = compute_rdm3s(self.det, self.coef[:, i], self.system.norbitals, self.system.noccupied_alpha, self.system.noccupied_beta)
            self.rdms[i]['aaa'] = dm3aaa
            self.rdms[i]['aab'] = dm3aab
            self.rdms[i]['abb'] = dm3abb
            self.rdms[i]['bbb'] = dm3bbb


