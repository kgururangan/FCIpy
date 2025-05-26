import numpy as np
from pyscf import gto, scf, fci
from fcipy.driver import Driver 

def test_fci_h4():

        Re = 1.0

        # As you make this more square, I think the diagonal initial guess will fail
        geom = [['H', (-Re, -Re, 0.000)],
                ['H', (-Re,  Re, 0.000)],
                ['H', (Re, -Re, 0.000)],
                ['H', (Re,  Re, 0.000)]]

        mol = gto.M(atom=geom, basis="dz", symmetry="D2H", spin=0, charge=0, unit="Bohr")
        mf = scf.RHF(mol)
        mf.kernel()

        #
        # create an FCI solver based on the SCF object
        #
        cisolver = fci.FCI(mf)
        cisolver.kernel(nroots=10)

        driver = Driver.from_pyscf(mf, nfrozen=0)
        driver.system.print_info()

        # obtain the determinant list for FCI
        driver.load_determinants(max_excit_rank=-1, target_irrep="AG")

        # perform a dense diagonalization of the full FCI Hamiltonian
        driver.build_hamiltonian(herm=True)
        driver.diagonalize_hamiltonian(herm=True)

        for i, e in enumerate(driver.total_energy[:20]):
               print(f"root {i}, E = {e}")
        print("FCI energies from PySCF:", cisolver.e_tot)

        #
        # Check the results
        #
        assert np.allclose(driver.total_energy[0], cisolver.e_tot[0], atol=1.0e-07)


if __name__ == "__main__":
        test_fci_h4()
