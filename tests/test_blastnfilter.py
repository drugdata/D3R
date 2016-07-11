#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_blastnfilter
----------------------------------

Tests for `blastnfilter` module.
"""

import unittest

from d3r import blastnfilter


class TestBlastnfilter(unittest.TestCase):

    def setUp(self):
        pass

    def test_parse_arguments(self):
        theargs = []
        try:
            blastnfilter._parse_arguments('hi', theargs)
            self.fail('expected exception')
        except:
            pass

        theargs = ['--outdir', 'myoutdir',
                   '--nonpolymertsv', 'mynonpolymertsv',
                   '--sequencetsv', 'mysequencetsv',
                   '--crystalpH', 'mycrystalph',
                   '--pdbblastdb', 'mypdbblastdb',
                   '--pdbdb', 'mypdbdb',
                   '--compinchi', 'mycompinchi']
        result = blastnfilter._parse_arguments('hi', theargs)
        self.assertEqual(result.out, 'myoutdir')
        self.assertEqual(result.non_polymer, 'mynonpolymertsv')
        self.assertEqual(result.loglevel, 'WARNING')
        self.assertEqual(result.polymer, 'mysequencetsv')
        self.assertEqual(result.ph, 'mycrystalph')
        self.assertEqual(result.blast_db, 'mypdbblastdb')
        self.assertEqual(result.pdb_path, 'mypdbdb')
        self.assertEqual(result.compinchi, 'mycompinchi')
        theargs.append('--log')
        theargs.append('DEBUG')
        result = blastnfilter._parse_arguments('hi', theargs)
        self.assertEqual(result.loglevel, 'DEBUG')

    def test_main_args_set_but_invalid(self):
        theargs = ['blastnfilter',
                   '--outdir', 'myoutdir',
                   '--nonpolymertsv', 'mynonpolymertsv',
                   '--sequencetsv', 'mysequencetsv',
                   '--crystalpH', 'mycrystalph',
                   '--pdbblastdb', 'mypdbblastdb',
                   '--pdbdb', 'mypdbdb',
                   '--compinchi', 'mycompinchi']
        try:
            blastnfilter.main(theargs)
            self.fail('expected exception')
        except:
            pass

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
