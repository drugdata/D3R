__author__ = 'churas'

import unittest
import tempfile
import os.path

from d3r.filter.filter import QueryFilter
"""
test_in_put
--------------------------------

Tests for `in_put` module.
"""

import shutil
from d3r.utilities.analysis import OutputAnalysis


class TestAnalysis(unittest.TestCase):
    def setUp(self):
        pass

    def test_print_filter_criteria(self):
        temp_dir = tempfile.mkdtemp()
        try:
            oa = OutputAnalysis()
            oa.print_filter_criteria(temp_dir)
            f = open(os.path.join(temp_dir, 'summary.txt'), 'r')
            lines = f.readlines()
            f.close()
            self.assertTrue('FILTERING CRITERIA' in lines[0])
            self.assertTrue('No. of query sequences           <=    ' +
                            str(QueryFilter.sequence_threshold) in lines[1])
            self.assertTrue('No. of dockable ligands           =    ' +
                            str(QueryFilter.dockable_ligand_threshold) in
                            lines[2])
            self.assertTrue('Ligand self symmetries' in lines[7])
            self.assertTrue('<=    100' in lines[7])
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
