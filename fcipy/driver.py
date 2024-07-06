"""Main calculation driver module of CIpy."""

import numpy as np
from importlib import import_module
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

    def load_determinants(self, file=None, method=None, target_irrep=None):
        from fcipy.determinants import read_determinants
        if file is not None:
            self.det = read_determinants(file)
        else:
            # import the specific CC method module and get its update function
            mod = import_module("fcipy.determinants")
            det_function = getattr(mod, method.lower())
            self.det = det_function(self.system, target_irrep)

        assert self.N_int == self.det.shape[0]
        self.ndet = self.det.shape[2]

    def run_ci(self, nroot):
        from fcipy.davidson import run_davidson
        self.coef = np.zeros((self.ndet, nroot))
        self.total_energy = np.zeros(nroot)
        self.total_energy, self.coef = run_davidson(self.system, self.det, self.e1int, self.e2int, nroot)

    def build_hamiltonian(self):
        from fcipy.lib import ci
        self.Hmat = ci.ci.build_hamiltonian(self.det, self.e1int, self.e2int,
                                            self.system.noccupied_alpha, self.system.noccupied_beta)

    def diagonalize_hamiltonian(self):
        if self.Hmat is None:
            self.build_hamiltonian()
        self.total_energy, self.coef = np.linalg.eigh(self.Hmat)
        self.total_energy += self.system.nuclear_repulsion



