#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'


import unittest

"""
test_makeblastdb
--------------------------------

Tests for `makeblastdb` module.
"""

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from tests.celpp import test_task


class TestMakeBlastDBTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_update_status_from_filesystem(self):
        params = D3RParameters()
        task = MakeBlastDBTask(None, params)
        test_task.try_update_status_from_filesystem(self, task)

    def test_constructor(self):
        params = D3RParameters()
        task = MakeBlastDBTask('/foo', params)
        self.assertEqual(task.get_name(), 'makeblastdb')
        self.assertEqual(task.get_stage(), 1)
        self.assertEqual(task.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(task.get_path(), '/foo')
        self.assertEqual(task.get_dir_name(), 'stage.1.makeblastdb')
        test_task.try_update_status_from_filesystem(self, task)

# Test run where pdbsequrl is not set

# Test run where makeblastdb is not set

# Test run where url download fails

# Test run where gunzip fails

# Test run where makeblastdb fails

# Test run where everythin is successful

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
