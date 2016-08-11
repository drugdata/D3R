__author__ = 'churas'

"""
test_vinadocking
--------------------------------

Tests for `vinadocking` module.
"""

import unittest
import os
import stat
import tempfile
import shutil

from d3r import chimera_proteinligprep


class fakeparams(object):
    """Holds fake parameters
    """
    pass


class TestChimera_ProteinLigPrep(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_ligand_prepare_ligfile_already_prepared(self):
        temp_dir = tempfile.mkdtemp()
        try:
            ligfile = os.path.join(temp_dir, 'ligfile')
            open(ligfile, 'a').close()
            val = chimera_proteinligprep.ligand_prepare(None, ligfile, None)
            self.assertTrue(val)
        finally:
            shutil.rmtree(temp_dir)

    def test_ligand_prepare_with_fake_rdkit_python(self):
        temp_dir = tempfile.mkdtemp()
        curdir = os.getcwd()
        try:
            os.chdir(temp_dir)
            ligfile = os.path.join(temp_dir, 'ligfile')
            lig_smile = os.path.join(temp_dir, '1fcz.smi')

            temp_bin_dir = os.path.join(temp_dir, 'bin')
            os.makedirs(temp_bin_dir)

            fakepython = os.path.join(temp_bin_dir, 'python')

            f = open(fakepython, 'w')
            f.write('#! /usr/bin/env python\n\n')
            f.write('import sys\n')
            f.write('f = open(sys.argv[2], "w")\n')
            f.write('f.write(sys.argv[1] + ":")\n')
            f.write('f.write(sys.argv[2] + ":")\n')
            f.write('f.write(sys.argv[3])\n')
            f.write('f.close()\n')
            f.write('sys.exit(0)\n')
            f.flush()
            f.close()
            os.chmod(fakepython, stat.S_IRWXU)

            val = chimera_proteinligprep.ligand_prepare(lig_smile,
                                                        ligfile, temp_dir)
            self.assertFalse(val)
            f = open(lig_smile)
            line = f.readline()
            unprep = lig_smile.replace('.smi','_unprep_step1.sdf')
            self.assertEqual(line, 'rdkit_smiles_to_3d_sdf.py:' +
                             lig_smile + ':' + unprep)
        finally:
            os.chdir(curdir)
            shutil.rmtree(temp_dir)

            
if __name__ == '__main__':
    unittest.main()
