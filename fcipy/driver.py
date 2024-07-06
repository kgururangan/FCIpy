"""Main calculation driver module of CIpy."""

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

    def run_ci(self, nroot, convergence=1.0e-08, max_size=30, maxit=200):
        from fcipy.davidson import run_davidson
        self.total_energy, self.coef = run_davidson(self.system, self.det, self.e1int, self.e2int, nroot,
                                                    convergence=convergence, max_size=max_size, maxit=maxit)

    def build_hamiltonian(self):
        from fcipy.hamiltonian import build_hamiltonian
        self.Hmat = build_hamiltonian(self.det, self.e1int, self.e2int, self.system.noccupied_alpha, self.system.noccupied_beta)

    def diagonalize_hamiltonian(self):
        if self.Hmat is None:
            self.build_hamiltonian()
        self.total_energy, self.coef = np.linalg.eigh(self.Hmat)
        self.total_energy += self.system.nuclear_repulsion



