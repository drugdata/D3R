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
            task = D3RTask(temp_dir, params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                            D3RTask.START_FILE)), True)
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.START_STATUS)
            task.start()
            self.assertNotEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask
                                                         .ERROR_FILE)), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_D3RTask_end(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(temp_dir, params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
            task.end()
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.
                                                         COMPLETE_FILE)), True)
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.COMPLETE_STATUS)
            task.set_error('some error')
            task.end()
            self.assertEqual(task.get_error(), 'some error')
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                            D3RTask.ERROR_FILE)), True)
            task.set_status(D3RTask.ERROR_STATUS)
            os.remove(os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE))
            task.set_error(None)
            task.end()
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
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
        self.assertEqual(task.get_nonpolymer_tsv(),
                         '/foo/stage.1.dataimport' +
                         '/new_release_structure_nonpolymer.tsv')

    def test_DataImportTask_get_sequence_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_sequence_tsv(),
                         '/foo/stage.1.dataimport/' +
                         'new_release_structure_sequence.tsv')

    def test_DataImportTask_get_crystalph_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_crystalph_tsv(),
                         '/foo/stage.1.dataimport/' +
                         'new_release_crystallization_pH.tsv')

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

    def test_BlastnFilterTask_get_txt_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            blast_task = BlastNFilterTask(temp_dir, params)
            blast_task.create_dir()
            # try empty directory
            self.assertEquals(len(blast_task.get_txt_files()), 0)

            # try where only summary.txt exists
            summary_file = os.path.join(blast_task.get_dir(), 'summary.txt')
            open(summary_file, 'a').close()
            self.assertEquals(len(blast_task.get_txt_files()), 0)

            # try non txt files only
            os.remove(summary_file)
            open(os.path.join(blast_task.get_dir(), 'foo.csv'), 'a').close()
            self.assertEquals(len(blast_task.get_txt_files()), 0)

            # try 1 txt file
            open(os.path.join(blast_task.get_dir(), '4vwx.txt'), 'a').close()
            self.assertEquals(len(blast_task.get_txt_files()), 1)

            # try 1 summary.txt and 1 txt file
            summary_file = os.path.join(blast_task.get_dir(), 'summary.txt')
            open(summary_file, 'a').close()
            self.assertEquals(len(blast_task.get_txt_files()), 1)

            # try multiple txt files
            open(os.path.join(blast_task.get_dir(), '5vwx.txt'), 'a').close()
            self.assertEquals(len(blast_task.get_txt_files()), 2)

        finally:
            shutil.rmtree(temp_dir)

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
            self.assertEqual(blastTask.get_error(), 'dataimport task has ' +
                             'notfound status')

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

            # try where compinchi is not complete
            completeFile = os.path.join(dataImport.get_path(),
                                        dataImport.get_dir_name(),
                                        D3RTask.COMPLETE_FILE)
            open(completeFile, 'a').close()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(), 'compinchi task has ' +
                             'notfound status')

            # try where compinchi failed
            compinchi = CompInchiDownloadTask(tempDir, params)
            compinchi.create_dir()
            errorFile = os.path.join(compinchi.get_path(),
                                     compinchi.get_dir_name(),
                                     D3RTask.ERROR_FILE)
            open(errorFile, 'a').close()
            self.assertEqual(blastTask.can_run(), False)
            self.assertEqual(blastTask.get_error(),
                             'compinchi task has error status')

            # try where blast can run
            completeFile = os.path.join(compinchi.get_path(),
                                        compinchi.get_dir_name(),
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
            params.blastdir = temp_dir
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
            f = open(std_out_file, 'r')
            echo_out = f.read().replace('\n', '')
            echo_out.index('--nonpolymertsv ' +
                           os.path.join(temp_dir, 'stage.1.dataimport',
                                        'new_release_structure_nonpolymer.tsv'
                                        ))
            echo_out.index('--sequencetsv ' +
                           os.path.join(temp_dir, 'stage.1.dataimport',
                                        'new_release_structure_sequence.tsv'))
            echo_out.index('--pdbblastdb ' +
                           os.path.join(temp_dir, 'current'))
            echo_out.index('--compinchi ' +
                           os.path.join(temp_dir, 'stage.1.compinchi',
                                        'Components-inchi.ich'))
            echo_out.index('--outdir ' +
                           os.path.join(temp_dir, 'stage.2.blastnfilter'))
            f.close()

            self.assertEqual(os.path.isfile(std_out_file), True)

            res = blasttask.get_email_log().rstrip('\n')
            res.index('/bin/echo')
            res.index('# targets found: 0')
            res.index('Target')
        finally:
            shutil.rmtree(temp_dir)

    def test_BlastNFilterTask_run_success_with_txt_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            foo_script = os.path.join(temp_dir, 'foo.py')
            params.blastnfilter = foo_script
            params.blastdir = temp_dir
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True

            txt_file = os.path.join(blasttask.get_dir(), '4za4.txt')

            txt_contents = ('query, 4za4\\n' +
                            'ph, 7.4\\n' +
                            'ligand, 4LU\\n\\n' +
                            'inchi, InChI=1S/C22H29N4O9P/c1-10-7-12-16-' +
                            '15(11(10)2)22(3,4)5-6-25(16)17-19(23-21(31)' +
                            '24-20(17)30)26(12)8-13(27)18(29)14(28)9-35-' +
                            '36(32,33)34/h6-7,13-14,18,27-29H,5,8-9H2,1-' +
                            '4H3,(H3-,23,24,30,31,32,33,34)/p+1/t13-,14+' +
                            ',18-/m0/s1\\n'
                            'largest, 4zz3, 4PP\\n' +
                            'smallest, 3ax3, 4LP\\n' +
                            'holo, 2ax1, XDN\\n' +
                            'apo, 2ll3, GSS\\n')

            # create fake blastnfilter script that makes csv files
            f = open(foo_script, 'w')
            f.write('#! /usr/bin/env python\n\n')
            f.write('f = open(\'' + txt_file + '\', \'w\')\n')
            f.write('f.write(\'' + txt_contents + '\\n\')\n')
            f.write('f.flush()\nf.close()\n')
            f.flush()
            f.close()
            os.chmod(foo_script, stat.S_IRWXU)

            blasttask.run()
            self.assertEqual(blasttask.get_error(), None)
            self.assertEqual(blasttask.get_status(), D3RTask.COMPLETE_STATUS)
            complete_file = os.path.join(blasttask.get_dir(),
                                         D3RTask.COMPLETE_FILE)

            self.assertEqual(os.path.isfile(complete_file), True)

            std_err_file = os.path.join(blasttask.get_dir(),
                                        'foo.py.stderr')

            self.assertEqual(os.path.isfile(std_err_file), True)

            std_out_file = os.path.join(blasttask.get_dir(),
                                        'foo.py.stdout')

            self.assertEqual(os.path.isfile(std_out_file), True)

            res = blasttask.get_email_log().rstrip('\n')
            res.index('/foo.py')
            res.index('# targets found: 1')
            res.index('Target')
            res.index('4za4')
        finally:
            shutil.rmtree(temp_dir)

    def test_BlastNFilterTask_run_with_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = 'false'
            params.blastdir = temp_dir
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
            params.blastdir = temp_dir
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
        params.blastdir = '/foo'
        blasttask = BlastNFilterTask(None, params)
        blasttask._can_run = False
        blasttask.run()

    def test_parse_blastnfilter_output_for_hit_stats(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.blastdir = temp_dir
            blasttask = BlastNFilterTask(temp_dir, params)

            blasttask.create_dir()
            self.assertEquals(blasttask
                              ._parse_blastnfilter_output_for_hit_stats(),
                              '\n# targets found: 0\n\n\nTarget\n')

            csv_file = os.path.join(blasttask.get_dir(), '4zyc.txt')
            f = open(csv_file, 'w')
            f.write('NEED TO PUT REAL DATA IN HERE\n')
            f.flush()
            f.close()
            self.assertEquals(blasttask
                              ._parse_blastnfilter_output_for_hit_stats(),
                              '\n# targets found: 1\n\n\nTarget\n' +
                              '4zyc\n')

            csv_file = os.path.join(blasttask.get_dir(), '4qqq.txt')
            f = open(csv_file, 'w')
            f.write('NEED TO PUT REAL DATA IN HERE\n')
            f.flush()
            f.close()
            res = blasttask._parse_blastnfilter_output_for_hit_stats()\
                .rstrip('\n')
            res.index('# targets found: 2')
            res.index('Target')
            res.index('4qqq')
            res.index('4zyc')

            csv_file = os.path.join(blasttask.get_dir(), '4abc.txt')
            f = open(csv_file, 'w')
            f.write('NEED TO PUT REAL DATA IN HERE\n')
            f.flush()
            f.close()
            res = blasttask._parse_blastnfilter_output_for_hit_stats()\
                .rstrip('\n')
            res.index('# targets found: 3')
            res.index('Target')
            res.index('4qqq')
            res.index('4zyc')
            res.index('4abc')

        finally:
            shutil.rmtree(temp_dir)

    def test_PDBPrepTask_can_run(self):
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

    def test_PDBPrepTask_run(self):
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

    def test_CompInchiDownloadTask_name_and_get_components_inchi_file(self):

        params = D3RParameters()
        task = CompInchiDownloadTask('foo', params)
        self.assertEquals(task.get_dir_name(), 'stage.1.compinchi')
        self.assertEquals(task.get_components_inchi_file(),
                          os.path.join('foo', 'stage.1.compinchi',
                                       'Components-inchi.ich'))

    def test_CompInchiDownloadTask_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # task does not exist
            params = D3RParameters()
            task = CompInchiDownloadTask(temp_dir, params)
            self.assertEquals(task.can_run(), True)

            task.create_dir()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(), task.get_dir_name() +
                              ' already exists and status is ' +
                              D3RTask.UNKNOWN_STATUS)

            open(os.path.join(task.get_dir(), D3RTask.START_FILE),
                 'a').close()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(), task.get_dir_name() +
                              ' already exists and status is ' +
                              D3RTask.START_STATUS)

            open(os.path.join(task.get_dir(), D3RTask.ERROR_FILE),
                 'a').close()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(), task.get_dir_name() +
                              ' already exists and status is ' +
                              D3RTask.ERROR_STATUS)

            os.remove(os.path.join(task.get_dir(), D3RTask.ERROR_FILE))

            open(os.path.join(task.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(), None)
        finally:
            shutil.rmtree(temp_dir)

    def test_CompInchiDownloadTask_run_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                # url not even set
                params = D3RParameters()
                task = CompInchiDownloadTask(temp_dir, params)
                task.run()
                self.fail('Expected exception')
            except AttributeError:
                pass

            shutil.rmtree(task.get_dir())
            # url invalid
            params = D3RParameters()
            params.compinchi = 'file:///' + os.path.join(temp_dir,
                                                         'doesnotexistasdf')
            task = CompInchiDownloadTask(temp_dir, params)
            task.run()
            self.assertEquals(task.get_error(), 'Unable to download file ' +
                              'from ' + params.compinchi)
        finally:
            shutil.rmtree(temp_dir)

    def test_CompInchiDownloadTask_run_success(self):
        temp_dir = tempfile.mkdtemp()
        try:

            params = D3RParameters()
            fakefile = os.path.join(temp_dir, 'fakefile')
            f = open(fakefile, 'w')
            f.write('hi\n')
            f.flush()
            f.close()

            params.compinchi = 'file://' + fakefile

            task = CompInchiDownloadTask(temp_dir, params)

            task.run()
            self.assertEquals(task.get_error(), None)
            os.path.isfile(task.get_components_inchi_file())
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
