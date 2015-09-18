#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import os.path
from d3r.filter.Filter import Filter
from d3r.Blast.Target import Target

"""
test_task
----------------------------------

Tests for `task` module.
"""

import shutil

from d3r import util



class TestUtil(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_filtered_targets(self):
        foo = Filter(None)
        self.assertEqual(foo.get_filtered_targets(), None)
        foo = Filter('cat')
        self.assertEqual(foo.get_filtered_targets(), 'cat')

    def test_filter_test_structures_by_identity(self):
        bar = Filter(None)
        bar.filter_test_structures_by_identity()


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
