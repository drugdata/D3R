#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_blastnfilter
----------------------------------

Tests for `blastnfilter` module.
"""

import unittest
import tempfile
import logging
import os
import os.path
import shutil
from datetime import date
from d3r import blastnfilter


class TestBlastnfilter(unittest.TestCase):

    def setUp(self):
        pass

    def test_setup_logging(self):
        logger = logging.getLogger('funlogger')

        blastnfilter._setup_logging(logging.INFO, blastnfilter.LOG_FORMAT)
        self.assertEqual(logging.getLogger('d3r.blastnfilter')
                         .getEffectiveLevel(),
                         logging.INFO)
        self.assertEqual(logging.getLogger('d3r.blast.ligand')
                         .getEffectiveLevel(),
                         logging.INFO)
        blastnfilter._setup_logging(logging.DEBUG, blastnfilter.LOG_FORMAT)
        self.assertEqual(logging.getLogger('d3r.blastnfilter')
                         .getEffectiveLevel(),
                         logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.blast.ligand')
                         .getEffectiveLevel(),
                         logging.DEBUG)


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
