__author__ = 'churas'

import unittest


"""
test_run
--------------------------------

Tests for `run` module.
"""

from d3r.utilities import run
from d3r import blastnfilter


class TestBlastNFilterTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_split_input(self):
        theargs = ['--outdir', '/f/myoutdir',
                   '--nonpolymertsv', '/f/mynonpolymertsv',
                   '--sequencetsv', '/f/mysequencetsv',
                   '--crystalpH', '/f/mycrystalph',
                   '--pdbblastdb', '/f/mypdbblastdb',
                   '--pdbdb', '/f/mypdbdb',
                   '--compinchi', '/f/mycompinchi']
        result = blastnfilter._parse_arguments('hi', theargs)
        non_polymer, polymer, ph, out_dir, blast_dir, pdb_db, pdb_path,\
            fasta, compinchi = run.split_input(result)
        self.assertEqual(non_polymer, '/f/mynonpolymertsv')
        self.assertEqual(polymer, '/f/mysequencetsv')
        self.assertEqual(ph, '/f/mycrystalph')
        self.assertEqual(out_dir, '/f/myoutdir')
        self.assertEqual(blast_dir, '/f/mypdbblastdb')
        self.assertEqual(pdb_db, '/f/mypdbblastdb/pdb_db')
        self.assertEqual(pdb_path, '/f/mypdbdb')
        self.assertEqual(fasta, '/f/mypdbblastdb/pdb_seqres.txt')
        self.assertEqual(compinchi, '/f/mycompinchi')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
