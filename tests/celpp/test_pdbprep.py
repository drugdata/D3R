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
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.pdbprep import PDBPrepTask


class TestPDBPrepTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no blast task found so it cannot run
            params = D3RParameters()
            pdbpreptask = PDBPrepTask(temp_dir, params)
            self.assertEqual(pdbpreptask.can_run(), False)
            self.assertEqual(pdbpreptask.get_error(),
                             'blastnfilter task has notfound status')

            # blastn filter running
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.START_FILE),
                 'a').close()
            pdbpreptask = PDBPrepTask(temp_dir, params)
            self.assertEqual(pdbpreptask.can_run(), False)
            self.assertEqual(pdbpreptask.get_error(),
                             'blastnfilter task has start status')

            # blastnfilter failed
            error_file = os.path.join(blastnfilter.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            pdbpreptask = PDBPrepTask(temp_dir, params)
            self.assertEqual(pdbpreptask.can_run(), False)
            self.assertEqual(pdbpreptask.get_error(),
                             'blastnfilter task has error status')

            # blastnfilter success
            os.remove(error_file)
            open(os.path.join(blastnfilter.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            pdbpreptask = PDBPrepTask(temp_dir, params)
            self.assertEqual(pdbpreptask.can_run(), True)
            self.assertEqual(pdbpreptask.get_error(), None)

            # pdbprep task exists already
            pdbpreptask = PDBPrepTask(temp_dir, params)
            pdbpreptask.create_dir()
            self.assertEqual(pdbpreptask.can_run(), False)
            self.assertEqual(pdbpreptask.get_error(),
                             pdbpreptask.get_dir_name() +
                             ' already exists and status is unknown')

            # pdbprep already complete
            pdbpreptask = PDBPrepTask(temp_dir, params)
            open(os.path.join(pdbpreptask.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(pdbpreptask.can_run(), False)
            self.assertEqual(pdbpreptask.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbprep = 'echo'

            # return immediately cause can_run is false
            pdbpreptask = PDBPrepTask(temp_dir, params)
            pdbpreptask.run()
            self.assertEqual(pdbpreptask.get_error(),
                             'blastnfilter task has notfound status')

            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
