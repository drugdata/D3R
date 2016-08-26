__author__ = 'churas'

"""
test_vinadocking
--------------------------------

Tests for `vinadocking` module.
"""

import unittest

from d3r import vinadocking


class TestVinaDocking(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_dock_center_missing_comma(self):
        center = '1,2'
        try:
            vinadocking.dock('ligand', 'protein', center)
        except IndexError as ie:
            self.assertEqual(str(ie), 'list index out of range', str(ie))

if __name__ == '__main__':
    unittest.main()
