#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import os.path

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

    def test_find_latest_year(self):
        tempDir = tempfile.mkdtemp()
        try:
            self.assertEqual(util.find_latest_year(tempDir), None)

            os.mkdir(os.path.join(tempDir, "foo"))
            os.mkdir(os.path.join(tempDir, "2014"))
            self.assertEqual(util.find_latest_year(tempDir),
                             os.path.join(tempDir, '2014'))

            open(os.path.join(tempDir, '2016'), 'a').close()
            self.assertEqual(util.find_latest_year(tempDir),
                             os.path.join(tempDir, '2014'))

            os.mkdir(os.path.join(tempDir, '2012'))
            self.assertEqual(util.find_latest_year(tempDir),
                             os.path.join(tempDir, '2014'))

            os.mkdir(os.path.join(tempDir, '2015'))

            self.assertEqual(util.find_latest_year(tempDir),
                             os.path.join(tempDir, '2015'))

        finally:
            shutil.rmtree(tempDir)

    def test_find_latest_weekly_dataset(self):
        tempDir = tempfile.mkdtemp()

        try:
            self.assertEqual(util.find_latest_weekly_dataset(tempDir), None)

            os.mkdir(os.path.join(tempDir, '2015'))

            self.assertEqual(util.find_latest_weekly_dataset(tempDir), None)

            open(os.path.join(tempDir, '2015', 'dataset.week.4'), 'a').close()

            self.assertEqual(util.find_latest_weekly_dataset(tempDir), None)

            os.mkdir(os.path.join(tempDir, '2015', 'dataset.week.3'))

            self.assertEqual(util.find_latest_weekly_dataset(tempDir),
                             os.path.join(tempDir, '2015',
                                          'dataset.week.3'))
        finally:
            shutil.rmtree(tempDir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
