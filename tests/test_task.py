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


from d3r import task
from d3r.task import D3RParameters
from d3r.task import D3RTask
from d3r.task import BlastNFilterTask

class TestD3rTask(unittest.TestCase):

    def setUp(self):
        pass

    def test_find_latest_year(self):
        tempDir = tempfile.mkdtemp()
        self.assertEqual(task.find_latest_year(tempDir), None)

        os.mkdir(os.path.join(tempDir, "foo"))
        os.mkdir(os.path.join(tempDir, "2014"))
        self.assertEqual(task.find_latest_year(tempDir),
                         os.path.join(tempDir, '2014'))

        open(os.path.join(tempDir, '2016'), 'a').close()
        self.assertEqual(task.find_latest_year(tempDir),
                         os.path.join(tempDir, '2014'))

        os.mkdir(os.path.join(tempDir, '2012'))
        self.assertEqual(task.find_latest_year(tempDir),
                         os.path.join(tempDir, '2014'))

        os.mkdir(os.path.join(tempDir, '2015'))

        self.assertEqual(task.find_latest_year(tempDir),
                         os.path.join(tempDir, '2015'))

        os.rmdir(os.path.join(tempDir, "foo"))
        os.rmdir(os.path.join(tempDir, "2012"))
        os.rmdir(os.path.join(tempDir, "2014"))
        os.rmdir(os.path.join(tempDir, "2015"))
        os.remove(os.path.join(tempDir, '2016'))
        os.rmdir(tempDir)

    def test_find_latest_weekly_dataset(self):
        tempDir = tempfile.mkdtemp()

        self.assertEqual(task.find_latest_weekly_dataset(tempDir), None)

        os.mkdir(os.path.join(tempDir, '2015'))

        self.assertEqual(task.find_latest_weekly_dataset(tempDir), None)

        open(os.path.join(tempDir, '2015', 'dataset.week.4'), 'a').close()

        self.assertEqual(task.find_latest_weekly_dataset(tempDir), None)

        os.mkdir(os.path.join(tempDir, '2015', 'dataset.week.3'))

        self.assertEqual(task.find_latest_weekly_dataset(tempDir),
                         os.path.join(tempDir, '2015',
                                      'dataset.week.3'))

        os.rmdir(os.path.join(tempDir, "2015", "dataset.week.3"))
        os.remove(os.path.join(tempDir, "2015", "dataset.week.4"))
        os.rmdir(os.path.join(tempDir, "2015"))
        os.rmdir(tempDir)

    def test_D3RTask(self):
        params = D3RParameters()
        task = D3RTask('/path',params)
        self.assertEqual(task.get_name(), None)
        self.assertEqual(task.get_path(),'/path')
        self.assertEqual(task.get_stage(), None)
        self.assertEqual(task.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(task.get_error(), None)

        task.set_name('foo')
        task.set_path('blah')
        task.set_stage(4)
        task.set_status(D3RTask.START_STATUS)
        task.set_error('error')

        self.assertEqual(task.get_name(), 'foo')
        self.assertEqual(task.get_path(), 'blah')
        self.assertEqual(task.get_stage(), 4)
        self.assertEqual(task.get_status(), D3RTask.START_STATUS)
        self.assertEqual(task.get_error(), 'error')

    def test_BlastNFilterTask(self):
        params = D3RParameters()
        blasttask = BlastNFilterTask('ha',params)
        self.assertEqual(blasttask.get_name(),'blastnfilter')
        self.assertEqual(blasttask.get_path(),'ha')
        self.assertEqual(blasttask.get_stage(),2)
        self.assertEqual(blasttask.get_status(),D3RTask.UNKNOWN_STATUS)
        self.assertEqual(blasttask.get_error(),None)


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
