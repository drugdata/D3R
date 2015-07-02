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
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                fakefile = os.path.join(temp_dir,'fakefile')
                open(fakefile,'a').close()
                util.find_latest_year(fakefile)
                self.fail('Expected exception')
            except Exception:
                pass
            self.assertEqual(util.find_latest_year(temp_dir), None)

            os.mkdir(os.path.join(temp_dir, "foo"))
            os.mkdir(os.path.join(temp_dir, "2014"))
            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2014'))

            open(os.path.join(temp_dir, '2016'), 'a').close()
            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2014'))

            os.mkdir(os.path.join(temp_dir, '2012'))
            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2014'))

            os.mkdir(os.path.join(temp_dir, '2015'))

            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2015'))

        finally:
            shutil.rmtree(temp_dir)

    def test_find_latest_weekly_dataset(self):
        temp_dir = tempfile.mkdtemp()

        try:
            self.assertEqual(util.find_latest_weekly_dataset(temp_dir), None)

            os.mkdir(os.path.join(temp_dir, '2015'))

            self.assertEqual(util.find_latest_weekly_dataset(temp_dir), None)

            open(os.path.join(temp_dir, '2015', 'dataset.week.4'), 'a').close()

            self.assertEqual(util.find_latest_weekly_dataset(temp_dir), None)

            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.3'))

            self.assertEqual(util.find_latest_weekly_dataset(temp_dir),
                             os.path.join(temp_dir, '2015',
                                          'dataset.week.3'))
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
