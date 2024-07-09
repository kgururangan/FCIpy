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

    def __init__(self, system, e1int, e2int):
        self.system = system
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
        print(f"   Number of unique alpha strings: {self.num_alpha}")
        print(f"   Number of unique beta strings: {self.num_beta}")

    def print_determinants(self):
        for idet in range(self.ndet):
            print(f"determinant {idet + 1}: alpha = {self.det[:, 0, idet]}  beta = {self.det[:, 1, idet]}")

    def run_ci(self, nroot, convergence=1.0e-08, max_size=30, maxit=200, opt=True):
        from fcipy.davidson import run_davidson, run_davidson_opt
        if opt:
            self.total_energy, self.coef = run_davidson_opt(self.system, self.det, self.num_alpha, self.num_beta, self.e1int, self.e2int, nroot,
                                                       convergence=convergence, max_size=max_size, maxit=maxit)
        else:
            self.total_energy, self.coef = run_davidson(self.system, self.det, self.e1int, self.e2int, nroot,
                                                        convergence=convergence, max_size=max_size, maxit=maxit)

    def build_hamiltonian(self, opt=True):
        from fcipy.hamiltonian import build_hamiltonian, build_hamiltonian_opt
        if opt:
            self.det, self.Hmat = build_hamiltonian_opt(self.det, self.num_alpha, self.num_beta, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta)
        else:
            self.Hmat = build_hamiltonian(self.det, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta)

    def diagonalize_hamiltonian(self, opt=True):
        if self.Hmat is None:
            self.build_hamiltonian(opt=opt)
        self.total_energy, self.coef = np.linalg.eigh(self.Hmat)
        self.total_energy += self.system.nuclear_repulsion

    def one_e_density_matrix(self):
        from fcipy.density import compute_rdm1
        for i in range(self.coef.shape[1]):
            self.rdm1[i] = compute_rdm1(self.det, self.coef[:, i], self.system.norbitals)
            self.nat_occ_num[i], self.nat_orb[i] = np.linalg.eig(self.rdm1[i])
        

