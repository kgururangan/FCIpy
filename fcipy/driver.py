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

    def load_determinants(self, file):
        from fcipy.determinants import read_determinants
        self.det = read_determinants(file)
        self.ndet = self.det.shape[2]

    def run_ci(self, nroot):
        from fcipy.davidson import run_davidson
        self.coef = np.zeros((self.ndet, nroot))
        self.total_energy = np.zeros(nroot)
        self.total_energy, self.coef = run_davidson(self.system, self.det, self.e1int, self.e2int, nroot)



