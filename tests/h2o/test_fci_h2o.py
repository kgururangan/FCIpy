import numpy as np
from fcipy.driver import Driver

def test_fci_h2o():
    from pyscf import gto, scf, fci

    geom = [["H", (0, 1.515263, -1.058898)],
            ["H", (0, -1.515263, -1.058898)],
            ["O", (0.0, 0.0, -0.0090)]]

    mol = gto.M(atom=geom, basis="sto-3g", spin=0, charge=0, unit="Bohr", symmetry="C2V")
    mf = scf.RHF(mol)
    mf.kernel()

    #
    # create an FCI solver based on the SCF object
    #
    mc = fci.FCI(mf)
    expected_energy = mc.kernel()[0]
    (dm1a, dm1b), (dm2aa, dm2ab, dm2bb), (dm3aaa, dm3aab, dm3abb, dm3bbb) = fci.direct_spin1.make_rdm123s(mc.ci, mol.nao, (5, 5), reorder=True)
    # [WARNING]: 2-RDMs and 3-RDMs from PySCF must be transposed from Chemist to Physics notation
    dm2aa = dm2aa.transpose(0, 2, 1, 3)
    dm2ab = dm2ab.transpose(0, 2, 1, 3)
    dm2bb = dm2bb.transpose(0, 2, 1, 3)
    dm3aaa = dm3aaa.transpose(0, 2, 4, 1, 3, 5)
    dm3aab = dm3aab.transpose(0, 2, 4, 1, 3, 5)
    dm3abb = dm3abb.transpose(0, 2, 4, 1, 3, 5)
    dm3bbb = dm3bbb.transpose(0, 2, 4, 1, 3, 5)

    expected_rdms = {'a': dm1a, 'b': dm1b, 'aa': dm2aa, 'ab': dm2ab, 'bb': dm2bb, 'aaa': dm3aaa, 'aab': dm3aab, 'abb': dm3abb, 'bbb': dm3bbb}

    driver = Driver.from_pyscf(mf, nfrozen=0)
    driver.system.print_info()

    # read the determinants
    driver.load_determinants(max_excit_rank=-1, target_irrep="A1")

    # run CI calculation
    driver.run_ci(nroot=1, prtol=0.01, herm=True)

    # compute the 1-RDM
    driver.build_rdm1s(i=0)
    # compute the 2-RDM
    driver.build_rdm2s(i=0)
    # compute the 3-RDM
    driver.build_rdm3s(i=0)

    # check the results
    assert np.allclose(driver.total_energy[0], expected_energy)
    for key, value in driver.rdms[0].items():
        assert np.allclose(np.linalg.norm(value.flatten() - expected_rdms[key].flatten()), 0.0, atol=1.0e-08, rtol=1.0e-08)


if __name__ == "__main__":
    test_fci_h2o()
