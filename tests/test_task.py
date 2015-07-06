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
import platform
import os
from d3r.task import D3RParameters
from d3r.task import UnsetPathError
from d3r.task import UnsetStageError
from d3r.task import UnsetNameError
from d3r.task import UnsetFileNameError
from d3r.task import D3RTask
from d3r.task import BlastNFilterTask
from d3r.task import DataImportTask
from d3r.task import MakeBlastDBTask


class TestD3rTask(unittest.TestCase):

    def setUp(self):
        pass

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
            self.fail('Expected UnsetStageError')
        except UnsetStageError:
            pass

        task.set_stage(1)
        try:
            task.get_dir_name()
            self.fail('Expected UnsetNameError')
        except UnsetNameError:
            pass

        task.set_name('foo')
        self.assertEqual(task.get_dir_name(), 'stage.1.foo')

    def test_D3RTask_get_dir(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.get_dir()
            self.fail('Expected UnsetPathError')
        except UnsetPathError:
            pass

        task.set_path('/blah')

        try:
            task.get_dir()
            self.fail('Expected UnsetStageError')
        except UnsetStageError:
            pass

        task.set_stage(1)
        try:
            task.get_dir()
            self.fail('Expected UnsetNameError')
        except UnsetNameError:
            pass

        task.set_name('foo')

        self.assertEqual(task.get_dir(), '/blah/stage.1.foo')
    
    def test_D3RTask_can_run(self):
        task = D3RTask(None, D3RParameters())
        self.assertEqual(task.can_run(), False)

    def test_D3RTask_run(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        self.assertEqual(task._can_run, None)
        task.run()
        self.assertEqual(task._can_run, False)
        task._can_run = True
  
        temp_dir = tempfile.mkdtemp()
        task.set_name('foo')
        task.set_stage(1)
        task.set_path(temp_dir)
        try:
           task.run()
        finally:
            shutil.rmtree(temp_dir)

    def test_D3RTask_write_to_file(self):
        tempDir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            try:
                task.write_to_file('hello', None)
                self.fail('Expected UnsetFileNameError')
            except UnsetFileNameError:
                pass
            try:
                task.write_to_file('hello', 'foo')
                self.fail('Expected UnsetPathError')
            except UnsetPathError:
                pass
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(tempDir)

            try:
                task.write_to_file('hello', 'foo')
                self.fail('Expected IOError')
            except IOError:
                pass
            task.create_dir()
            task.write_to_file('hello', 'foo')
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                            'foo')), True)
        finally:
            shutil.rmtree(tempDir)

    def test_D3RTask_create_dir(self):
        tempDir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            try:
                task.create_dir()
                self.fail('Expected UnsetPathError')
            except UnsetPathError:
                pass
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
            self.fail("Expected UnsetPathError")
        except UnsetPathError:
            pass

        tempDir = tempfile.mkdtemp()
        try:
            task.set_path(tempDir)

            # Test unset stage
            try:
                task.update_status_from_filesystem()
                self.fail("Expected UnsetStageError")
            except UnsetStageError:
                pass
        finally:
            shutil.rmtree(tempDir)

        task.set_stage(1)
        task.set_name('foo')

        self._try_update_status_from_filesystem(task)

    def test_D3RTask_start(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(temp_dir,params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
	    self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.START_FILE)), True)
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(),D3RTask.START_STATUS) 
            task.start()
            self.assertNotEqual(task.get_error(), None)
            self.assertEqual(task.get_status(),D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.ERROR_FILE)), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_D3RTask_end(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(temp_dir,params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
            task.end()
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.COMPLETE_FILE)), True)
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(),D3RTask.COMPLETE_STATUS)
            task.set_error('some error')
            task.end()
            self.assertEqual(task.get_error(), 'some error')
            self.assertEqual(task.get_status(),D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.ERROR_FILE)), True)
            task.set_status(D3RTask.ERROR_STATUS)
            os.remove(os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE))
            task.set_error(None)
            task.end()
            self.assertEqual(task.get_status(),D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.ERROR_FILE)), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_D3RTask_get_smtp_server(self):
         params = D3RParameters()
         params.smtp = ''
         params.smtpport = '25'
         task = D3RTask(None, params)
         self.assertNotEqual(task._get_smtp_server(), None)

    def test_D3RTask_build_from_address(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        exp_from_addr = os.getlogin() + '@' + platform.node()
        self.assertEqual(task._build_from_address(), exp_from_addr)


    def test_BlastNFilterTask_update_status_from_filesystem(self):
        params = D3RParameters()
        task = BlastNFilterTask(None, params)
        self._try_update_status_from_filesystem(task)

    def test_DataImportTask_update_status_from_filesystem(self):
        params = D3RParameters()
        task = DataImportTask(None, params)
        self._try_update_status_from_filesystem(task)

    def test_DataImportTask_get_nonpolymer_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_nonpolymer_tsv(), '/foo/stage.1.dataimport/new_release_structure_nonpolymer.tsv')
        
    def test_DataImportTask_get_sequence_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_sequence_tsv(), '/foo/stage.1.dataimport/new_release_structure_sequence.tsv')
   
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

    def test_BlastNFilterTask_run_with_success(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = '/bin/echo'
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True
            blasttask.run()
            self.assertEqual(blasttask.get_status(), D3RTask.COMPLETE_STATUS)
            self.assertEqual(blasttask.get_error(), None)
            complete_file = os.path.join(blasttask.get_dir(),
                                         D3RTask.COMPLETE_FILE)

            self.assertEqual(os.path.isfile(complete_file), True)

            std_err_file = os.path.join(blasttask.get_dir(),
                                        'echo.stderr')

            self.assertEqual(os.path.isfile(std_err_file), True)

            std_out_file = os.path.join(blasttask.get_dir(),
                                        'echo.stdout')

            self.assertEqual(os.path.isfile(std_out_file), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_BlastNFilterTask_run_with_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = 'false'
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True
            blasttask.run()
            self.assertEqual(blasttask.get_status(), D3RTask.ERROR_STATUS)
            self.assertNotEqual(blasttask.get_error(), None)
            error_file = os.path.join(blasttask.get_dir(),
                                         D3RTask.ERROR_FILE)

            self.assertEqual(os.path.isfile(error_file), True)
            
            std_err_file = os.path.join(blasttask.get_dir(),
                                        'false.stderr')
            
            self.assertEqual(os.path.isfile(std_err_file), True)

            std_out_file = os.path.join(blasttask.get_dir(),
                                        'false.stdout')

            self.assertEqual(os.path.isfile(std_out_file), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_BlastNFilterTask_run_with_exception(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = 'falseasdfasdf'
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True
            blasttask.run()
            self.assertEqual(blasttask.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(blasttask.get_error().startswith('Caught'), True)
            self.assertNotEqual(blasttask.get_error(), None)
        finally:
            shutil.rmtree(temp_dir)


    def test_BlastNFilterTask_run_with_can_run_already_set_false(self):
        params = D3RParameters()
        params.blastnfilter = 'false'
        blasttask = BlastNFilterTask(None, params)
        blasttask._can_run = False
        blasttask.run()
   

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
