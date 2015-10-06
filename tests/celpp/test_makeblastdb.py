#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'


import unittest
import tempfile
import os.path

"""
test_blastnfilter
--------------------------------

Tests for `blastnfilter` module.
"""

import shutil
import platform
import os
import stat
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import UnsetPathError
from d3r.celpp.task import UnsetStageError
from d3r.celpp.task import UnsetNameError
from d3r.celpp.task import UnsetFileNameError
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp.pdbprep import PDBPrepTask
from d3r.celpp.compinchidownload import CompInchiDownloadTask

from tests.celpp import test_task

class TestMakeBlastDBTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_MakeBlastDBTask(self):
        params = D3RParameters()
        try:
            task = MakeBlastDBTask('blah', params)
            self.fail("Expected AttributeError")
        except AttributeError:
            pass
        params.blastdir = '/foo'
        task = MakeBlastDBTask('blah', params)
        self.assertEqual(task.get_name(), 'makeblastdb')
        self.assertEqual(task.get_stage(), 1)
        self.assertEqual(task.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(task.get_path(), '/foo')
        self.assertEqual(task.get_dir_name(), 'current')
        test_task.try_update_status_from_filesystem(self, task)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()