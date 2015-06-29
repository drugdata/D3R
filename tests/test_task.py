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

from d3r import task
from d3r.task import D3RParameters
from d3r.task import UnsetPathException
from d3r.task import UnsetStageException
from d3r.task import UnsetNameException
from d3r.task import D3RTask
from d3r.task import BlastNFilterTask
from d3r.task import DataImportTask
from d3r.task import MakeBlastDBTask


class TestD3rTask(unittest.TestCase):

    def setUp(self):
        pass

    def test_find_latest_year(self):
        tempDir = tempfile.mkdtemp()
        try:
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

        finally:
            shutil.rmtree(tempDir)

    def test_find_latest_weekly_dataset(self):
        tempDir = tempfile.mkdtemp()

        try:
            self.assertEqual(task.find_latest_weekly_dataset(tempDir), None)

            os.mkdir(os.path.join(tempDir, '2015'))

            self.assertEqual(task.find_latest_weekly_dataset(tempDir), None)

            open(os.path.join(tempDir, '2015', 'dataset.week.4'), 'a').close()

            self.assertEqual(task.find_latest_weekly_dataset(tempDir), None)

            os.mkdir(os.path.join(tempDir, '2015', 'dataset.week.3'))

            self.assertEqual(task.find_latest_weekly_dataset(tempDir),
                             os.path.join(tempDir, '2015',
                                          'dataset.week.3'))
        finally:
            shutil.rmtree(tempDir)

    def test_D3RTask(self):
        params = D3RParameters()
        task = D3RTask('/path', params)
        self.assertEqual(task.get_name(), None)
        self.assertEqual(task.get_path(), '/path')
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

    def test_get_dir_name(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.get_dir_name()
            self.fail('Expected UnsetStageException')
        except UnsetStageException:
            pass

        task.set_stage(1)
        try:
            task.get_dir_name()
            self.fail('Expected UnsetNameException')
        except UnsetNameException:
            pass

        task.set_name('foo')
        self.assertEqual(task.get_dir_name(), 'stage.1.foo')

    def test_D3RTask_create_dir(self):
        tempDir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            self.assertEqual(task.create_dir(), None)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(tempDir)
            self.assertEqual(task.create_dir(),
                             os.path.join(tempDir, 'stage.1.foo'))
        finally:
            shutil.rmtree(tempDir)

    def test_D3RTask_update_status_from_filesystem(self):
        params = D3RParameters()
        task = D3RTask(None, params)

        try:
            task.update_status_from_filesystem()
            self.fail("Expected UnsetPathException")
        except UnsetPathException:
            pass

        tempDir = tempfile.mkdtemp()
        try:
            task.set_path(tempDir)

            # Test unset stage
            try:
                task.update_status_from_filesystem()
                self.fail("Expected UnsetStageException")
            except UnsetStageException:
                pass
        finally:
            shutil.rmtree(tempDir)

        task.set_stage(1)
        task.set_name('foo')

        self._try_update_status_from_filesystem(task)

    def test_BlastNFilterTask_update_status_from_filesystem(self):
        params = D3RParameters()
        task = BlastNFilterTask(None, params)
        self._try_update_status_from_filesystem(task)

    def test_DataImportTask_update_status_from_filesystem(self):
        params = D3RParameters()
        task = DataImportTask(None, params)
        self._try_update_status_from_filesystem(task)

    def _try_update_status_from_filesystem(self, task):

        tempDir = tempfile.mkdtemp()
        try:
            task.set_path(tempDir)

            self.assertEqual(task.update_status_from_filesystem(),
                             D3RTask.NOTFOUND_STATUS)

            self.assertEqual(task.get_status(),
                             D3RTask.NOTFOUND_STATUS)

            # Test directory exists but no complete, start, or error files
            task.create_dir()
            self.assertEqual(task.update_status_from_filesystem(),
                             D3RTask.UNKNOWN_STATUS)

            self.assertEqual(task.get_status(),
                             D3RTask.UNKNOWN_STATUS)

            # Test start file exists
            startFile = os.path.join(tempDir, task.get_dir_name(),
                                     'start')
            open(startFile, 'a').close()

            self.assertEqual(task.update_status_from_filesystem(),
                             D3RTask.START_STATUS)

            self.assertEqual(task.get_status(),
                             D3RTask.START_STATUS)

            # Test error file exists
            errorFile = os.path.join(tempDir, task.get_dir_name(),
                                     'error')
            open(errorFile, 'a').close()

            self.assertEqual(task.update_status_from_filesystem(),
                             D3RTask.ERROR_STATUS)

            self.assertEqual(task.get_status(),
                             D3RTask.ERROR_STATUS)

            # Test complete file exists
            completeFile = os.path.join(tempDir, task.get_dir_name(),
                                        'complete')
            open(completeFile, 'a').close()

            self.assertEqual(task.update_status_from_filesystem(),
                             D3RTask.COMPLETE_STATUS)

            self.assertEqual(task.get_status(),
                             D3RTask.COMPLETE_STATUS)

        finally:
            shutil.rmtree(tempDir)

    def test_BlastNFilterTask(self):
        params = D3RParameters()
        blasttask = BlastNFilterTask('ha', params)
        self.assertEqual(blasttask.get_name(), 'blastnfilter')
        self.assertEqual(blasttask.get_path(), 'ha')
        self.assertEqual(blasttask.get_stage(), 2)
        self.assertEqual(blasttask.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(blasttask.get_error(), None)

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
        self._try_update_status_from_filesystem(task)

    def test_BlastNFilterTask_can_run(self):
        tempDir = tempfile.mkdtemp()

        try:
            # try where makeblastdb is not complete
            params = D3RParameters()
            params.blastdir = tempDir
            blastTask = BlastNFilterTask(tempDir, params)
            self.assertEqual(blastTask.can_run(), False)

            # try where makeblastdb failed
            blastDb = MakeBlastDBTask(tempDir, params)
            blastDb.create_dir()
            errorFile = os.path.join(blastDb.get_path(),
                                     blastDb.get_dir_name(),
                                     D3RTask.ERROR_FILE)
            open(errorFile, 'a').close()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(),
                             'makeblastdb task has error status')

            # try where data import is not complete
            completeFile = os.path.join(blastDb.get_path(),
                                        blastDb.get_dir_name(),
                                        D3RTask.COMPLETE_FILE)
            open(completeFile, 'a').close()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(), None)

            # try where data import failed
            dataImport = DataImportTask(tempDir, params)
            dataImport.create_dir()
            errorFile = os.path.join(dataImport.get_path(),
                                     dataImport.get_dir_name(),
                                     D3RTask.ERROR_FILE)
            open(errorFile, 'a').close()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(),
                             'dataimport task has error status')

            # try where blast can run
            completeFile = os.path.join(dataImport.get_path(),
                                        dataImport.get_dir_name(),
                                        D3RTask.COMPLETE_FILE)
            open(completeFile, 'a').close()
            self.assertEqual(blastTask.can_run(), True)
            self.assertEqual(blastTask.get_error(), None)

            # try where blast exists
            blastTask.create_dir()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(),
                             'stage.2.blastnfilter already exists and' +
                             ' status is unknown')

            # try where blast is complete
            completeFile = os.path.join(blastTask.get_path(),
                                        blastTask.get_dir_name(),
                                        D3RTask.COMPLETE_FILE)
            open(completeFile, 'a').close()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(), None)

        finally:
            shutil.rmtree(tempDir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
