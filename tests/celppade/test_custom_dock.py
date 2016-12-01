#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_custom_dock
----------------------------------

Tests for `custom_dock` module.
"""

import unittest
import tempfile
import shutil
import os

from d3r.celppade.custom_dock import Dock


class TestDock(unittest.TestCase):

    def setUp(self):
        pass

    def test_dock(self):
        """Test where we run dock() with fake arguments
           to verify we always get False
        """
        d = Dock()
        self.assertEqual(d.dock(1, 2, 3, 4), False)

    def test_ligand_technical_prep(self):
        """Test where we run ligand_technical_prep() with fake argument
           to verify we always get back the argument as an array
        """
        d = Dock()
        self.assertEqual(d.ligand_technical_prep('hi'), ['hi'])

    def test_receptor_technical_prep(self):
        """Test where we run receptor_technical_prep() with fake arguments
           to verify we always get back the 1st argument as an array
        """
        d = Dock()
        self.assertEqual(d.receptor_technical_prep('hi', 'foo'), ['hi'])

    def test_get_pocket_center(self):
        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            # test where file does not exist
            self.assertEqual(d.get_pocket_center(temp_dir), False)

            # test with valid center.txt file
            pcfile = os.path.join(temp_dir, 'center.txt')
            f = open(pcfile, 'w')
            f.write('1,2')
            f.flush()
            f.close()
            self.assertEqual(d.get_pocket_center(temp_dir), [1.0, 2.0])

            # test with INVALID center.txt file
            # TODO Fix get_pocket_center() to not raise error when contents of
            # TODO center.txt file are invalid
            # f = open(pcfile, 'w')
            # f.write('a')
            # f.flush()
            # f.close()
            # self.assertEqual(d.get_pocket_center(temp_dir), False)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
