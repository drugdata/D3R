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
import stat
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from tests.celpp import test_task


class TestBlastNFilterTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_update_status_from_filesystem(self):
        params = D3RParameters()
        task = BlastNFilterTask(None, params)
        test_task.try_update_status_from_filesystem(self, task)

    def test_constructor(self):
        params = D3RParameters()
        blasttask = BlastNFilterTask('ha', params)
        self.assertEqual(blasttask.get_name(), 'blastnfilter')
        self.assertEqual(blasttask.get_path(), 'ha')
        self.assertEqual(blasttask.get_stage(), 2)
        self.assertEqual(blasttask.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(blasttask.get_error(), None)

    def test_get_txt_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            blast_task = BlastNFilterTask(temp_dir, params)

            # try missing directory
            self.assertEquals(len(blast_task.get_txt_files()), 0)

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

    def test_get_uploadable_files(self):
        self.assertEqual('Need to ', 'test this')

    def test_can_run(self):
        tempDir = tempfile.mkdtemp()

        try:
            # try where makeblastdb is not complete
            params = D3RParameters()
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

            # try where blast can run
            os.remove(errorFile)
            completeFile = os.path.join(dataImport.get_dir(),
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

    def test_run_with_blast_success_postanalysis_fail(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = '/bin/echo'
            params.postanalysis = os.path.join(temp_dir, 'foo.py')
            params.pdbdb = '/pdbdb'
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True

            txt_file = os.path.join(blasttask.get_dir(), 'summary.txt')

            txt_contents = ('INPUT SUMMARY\\n' +
                            '  sequences:  177\\n' +
                            '  complexes:  149\\n')
            # create fake blastnfilter script that makes csv files
            f = open(params.postanalysis, 'w')
            f.write('#! /usr/bin/env python\n\n')
            f.write('f = open(\'' + txt_file + '\', \'w\')\n')
            f.write('f.write(\'' + txt_contents + '\\n\')\n')
            f.write('f.flush()\nf.close()\n')
            f.flush()
            f.close()
            os.chmod(params.postanalysis, stat.S_IRWXU)

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
            echo_out.index(' --sequencetsv ' +
                           os.path.join(temp_dir, 'stage.1.dataimport',
                                        'new_release_structure_sequence.tsv'))
            echo_out.index(' --pdbblastdb ' +
                           os.path.join(temp_dir, 'stage.1.makeblastdb'))
            echo_out.index(' --compinchi ' +
                           os.path.join(temp_dir, 'stage.1.dataimport',
                                        'Components-inchi.ich'))
            echo_out.index(' --outdir ' +
                           os.path.join(temp_dir, 'stage.2.blastnfilter'))
            echo_out.index(' --crystalpH ' +
                           os.path.join(temp_dir, 'stage.1.dataimport',
                                        'new_release_crystallization_pH.tsv'))
            echo_out.index(' --pdbdb /pdbdb ')
            f.close()

            self.assertEqual(os.path.isfile(std_out_file), True)
            self.assertEquals(blasttask.get_status(), D3RTask.COMPLETE_STATUS)
            self.assertEquals(os.path.exists(os.path.join(blasttask.get_dir(),
                                                          'foo.py.stderr')),
                              True)
            self.assertEquals(os.path.exists(os.path.join(blasttask.get_dir(),
                                                          'foo.py.stdout')),
                              True)
            res = blasttask.get_email_log().rstrip('\n')
            res.index('/bin/echo')
            res.index('# txt files found: 0')
            res.index('Output from summary.txt')
            res.index('  sequences:  177')
            res.index('  complexes:  149')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_with_blast_success_postanalysis_success_no_summary_file(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = 'true'
            params.postanalysis = '/bin/echo'
            params.pdbdb = '/pdbdb'
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
            echo_out.index('--compinchi ' +
                           os.path.join(temp_dir, 'stage.1.dataimport',
                                        'Components-inchi.ich'))
            echo_out.index(' ' + os.path.join(temp_dir,
                                              'stage.2.blastnfilter'))
            f.close()

            self.assertEqual(os.path.isfile(std_out_file), True)
            self.assertEquals(blasttask.get_status(), D3RTask.COMPLETE_STATUS)
            self.assertEquals(os.path.exists(os.path.join(blasttask.get_dir(),
                                                          'true.stderr')),
                              True)
            self.assertEquals(os.path.exists(os.path.join(blasttask.get_dir(),
                                                          'true.stdout')),
                              True)
            res = blasttask.get_email_log().rstrip('\n')
            res.index('/bin/echo')
            res.index('# txt files found: 0')
            res.index('Output from summary.txt')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_success_with_txt_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            foo_script = os.path.join(temp_dir, 'foo.py')
            params.blastnfilter = foo_script
            params.postanalysis = '/bin/echo'
            params.pdbdb = '/pdbdb'
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
            res.index('# txt files found: 1')
            res.index('Output from summary.txt')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_with_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = 'false'
            params.postanalysis = 'true'
            params.pdbdb = '/pdbdb'
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True
            blasttask.run()
            self.assertEqual(blasttask.get_status(), D3RTask.COMPLETE_STATUS)
            self.assertEqual(blasttask.get_error(), None)
            complete_file = os.path.join(blasttask.get_dir(),
                                         D3RTask.COMPLETE_FILE)

            self.assertEqual(os.path.isfile(complete_file), True)
            error_file = os.path.join(blasttask.get_dir(),
                                      D3RTask.ERROR_FILE)

            self.assertEqual(os.path.isfile(error_file), False)

            std_err_file = os.path.join(blasttask.get_dir(),
                                        'false.stderr')

            self.assertEqual(os.path.isfile(std_err_file), True)

            std_out_file = os.path.join(blasttask.get_dir(),
                                        'false.stdout')

            self.assertEqual(os.path.isfile(std_out_file), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_with_exception(self):
        temp_dir = tempfile.mkdtemp()

        try:
            params = D3RParameters()
            params.blastnfilter = 'falseasdfasdf'
            params.postanalysis = 'true'
            params.pdbdb = '/pdbdb'
            blasttask = BlastNFilterTask(temp_dir, params)
            blasttask._can_run = True
            blasttask.run()
            self.assertEqual(blasttask.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(blasttask.get_error().startswith('Caught'), True)
            self.assertNotEqual(blasttask.get_error(), None)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_with_can_run_already_set_false(self):
        params = D3RParameters()
        params.blastnfilter = 'false'
        params.postanalysis = 'false'
        params.pdbdb = '/pdbdb'
        blasttask = BlastNFilterTask(None, params)
        blasttask._can_run = False
        blasttask.run()

    def test_parse_blastnfilter_output_for_hit_stats(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            blasttask = BlastNFilterTask(temp_dir, params)

            blasttask.create_dir()

            # no summary.txt file
            self.assertEquals(blasttask
                              ._parse_blastnfilter_output_for_hit_stats(),
                              '\n# txt files found: 0\n\nOutput from ' +
                              'summary.txt\n')

            csv_file = os.path.join(blasttask.get_dir(), '4zyc.txt')
            f = open(csv_file, 'w')
            f.write('NEED TO PUT REAL DATA IN HERE\n')
            f.flush()
            f.close()
            self.assertEquals(blasttask
                              ._parse_blastnfilter_output_for_hit_stats(),
                              '\n# txt files found: 1\n\nOutput from ' +
                              'summary.txt\n')

            csv_file = os.path.join(blasttask.get_dir(), '4qqq.txt')
            f = open(csv_file, 'w')
            f.write('NEED TO PUT REAL DATA IN HERE\n')
            f.flush()
            f.close()
            res = blasttask._parse_blastnfilter_output_for_hit_stats()\
                .rstrip('\n')
            res.index('# txt files found: 2')
            res.index('Output from summary.txt')

            csv_file = os.path.join(blasttask.get_dir(), '4abc.txt')
            f = open(csv_file, 'w')
            f.write('NEED TO PUT REAL DATA IN HERE\n')
            f.flush()
            f.close()
            res = blasttask._parse_blastnfilter_output_for_hit_stats()\
                .rstrip('\n')
            res.index('# txt files found: 3')
            res.index('Output from summary.txt')

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
