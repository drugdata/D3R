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
from d3r.celpp.blastnfilter import BlastNFilterSummary
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
        self.assertEqual(blasttask.get_stage(), 3)
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
        temp_dir = tempfile.mkdtemp()
        try:
            # test where directory doesn't even exist
            params = D3RParameters()
            task = BlastNFilterTask(temp_dir, params)
            self.assertEqual(task.get_uploadable_files(), [])

            # test on empty dir
            task.create_dir()

            # test with blastn filter log
            logfile = os.path.join(task.get_dir(),
                                   BlastNFilterTask.BLASTNFILTER_LOG)
            open(logfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            self.assertEqual(flist[0], logfile)

            # test with dockable.xslx
            dockable = os.path.join(task.get_dir(),
                                    BlastNFilterTask.DOCKABLE_XSLX)
            open(dockable, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(dockable)

            # test with summary.txt
            summary = os.path.join(task.get_dir(),
                                   BlastNFilterTask.SUMMARY_TXT)
            open(summary, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(summary)

            # test with 1 txt file
            txt_one = os.path.join(task.get_dir(),
                                   '1fcz.txt')
            open(txt_one, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 4)
            flist.index(txt_one)

            # test with 2 txt files
            txt_two = os.path.join(task.get_dir(),
                                   '4asd.txt')
            open(txt_two, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 5)
            flist.index(txt_two)

            # test with additional stderr/stdout files
            errfile = os.path.join(task.get_dir(),
                                   'blastnfilter.py.stderr')
            open(errfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 6)
            flist.index(errfile)

            outfile = os.path.join(task.get_dir(),
                                   'blastnfilter.py.stdout')
            open(outfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 7)
            flist.index(outfile)
            flist.index(errfile)
            flist.index(txt_two)
            flist.index(txt_one)
            flist.index(summary)
            flist.index(dockable)
            flist.index(logfile)

        finally:
            shutil.rmtree(temp_dir)

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
                             blastTask.get_dir_name() + ' already exists and' +
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

    def test_run_with_blast_success_usenewsequencetsv(self):
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

            dtask = DataImportTask(temp_dir, params)
            dtask.create_dir()
            open(dtask.get_sequence_tsv(), 'a').close()

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

            dataimport = DataImportTask(temp_dir, params)
            makeblast = MakeBlastDBTask(temp_dir, params)

            f = open(std_out_file, 'r')
            echo_out = f.read().replace('\n', '')
            echo_out.index('--nonpolymertsv ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.NONPOLYMER_TSV
                                        ))
            echo_out.index(' --sequencetsv ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.SEQUENCE_TSV))
            echo_out.index(' --pdbblastdb ' +
                           os.path.join(temp_dir, makeblast.get_dir_name()))
            echo_out.index(' --compinchi ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.COMPINCHI_ICH))
            echo_out.index(' --outdir ' +
                           os.path.join(temp_dir, blasttask.get_dir_name()))
            echo_out.index(' --crystalpH ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.CRYSTALPH_TSV))
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
            self.assertEqual(res.find(dataimport.get_sequence_tsv() +
                                      ' file not found falling back to ' +
                                      dataimport.get_oldsequence_tsv()), -1)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_with_blast_success_useoldseq_and_postanalysis_fail(self):
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

            dataimport = DataImportTask(temp_dir, params)
            makeblast = MakeBlastDBTask(temp_dir, params)

            f = open(std_out_file, 'r')
            echo_out = f.read().replace('\n', '')
            echo_out.index('--nonpolymertsv ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.NONPOLYMER_TSV
                                        ))
            echo_out.index(' --sequencetsv ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.OLDSEQUENCE_TSV))
            echo_out.index(' --pdbblastdb ' +
                           os.path.join(temp_dir, makeblast.get_dir_name()))
            echo_out.index(' --compinchi ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.COMPINCHI_ICH))
            echo_out.index(' --outdir ' +
                           os.path.join(temp_dir, blasttask.get_dir_name()))
            echo_out.index(' --crystalpH ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.CRYSTALPH_TSV))
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
            res.index(dataimport.get_sequence_tsv() +
                      ' file not found falling back to ' +
                      dataimport.get_oldsequence_tsv())
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
            dataimport = DataImportTask(temp_dir, params)
            f = open(std_out_file, 'r')
            echo_out = f.read().replace('\n', '')
            echo_out.index('--compinchi ' +
                           os.path.join(temp_dir, dataimport.get_dir_name(),
                                        DataImportTask.COMPINCHI_ICH))

            echo_out.index(' ' + os.path.join(temp_dir,
                                              blasttask.get_dir_name()))
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

    def test_get_blastnfilter_summary(self):
        """Tests BlastNFilterTask.get_blastnfilter_summary call
        """
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            blasttask = BlastNFilterTask(temp_dir, params)
            summary = blasttask.get_blastnfilter_summary()
            self.assertEqual(summary.get_complexes(), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_blastnfilter_summary_none_for_path(self):
        summary = BlastNFilterSummary(None)
        self.assertEqual(summary.get_complexes(), 0)
        self.assertEqual(summary.get_dockable_complexes(), 0)
        self.assertEqual(summary.get_dockable_monomers(), 0)
        self.assertEqual(summary.get_targets_found(), 0)
        self.assertEqual(summary.get_week_number(), '0')
        self.assertEqual(summary.get_year(), '0')
        self.assertEqual(summary.get_csv(), '0,0,0,0,0,0')

    def test_blastnfilter_summary_no_summary_file(self):
        summary = BlastNFilterSummary('/foo')
        self.assertEqual(summary.get_complexes(), 0)
        self.assertEqual(summary.get_dockable_complexes(), 0)
        self.assertEqual(summary.get_dockable_monomers(), 0)
        self.assertEqual(summary.get_targets_found(), 0)
        self.assertEqual(summary.get_week_number(), '0')
        self.assertEqual(summary.get_year(), '0')
        self.assertEqual(summary.get_csv(), '0,0,0,0,0,0')

    def test_blastnfilter_summary_getters_and_setters(self):
        summary = BlastNFilterSummary('/foo')

        summary.set_complexes(5)
        summary.set_dockable_complexes(6)
        summary.set_dockable_monomers(7)
        summary.set_targets_found(8)

        self.assertEqual(summary.get_complexes(), 5)

        self.assertEqual(summary.get_dockable_complexes(), 6)
        self.assertEqual(summary.get_dockable_monomers(), 7)
        self.assertEqual(summary.get_targets_found(), 8)
        self.assertEqual(summary.get_csv(), '0,0,5,6,7,8')

    def test_blastnfilter_summary_week_and_year(self):

        blast = BlastNFilterTask('/foo', D3RParameters())
        summary = BlastNFilterSummary('/foo/2018/dataset.week.4'
                                      '/' + blast.get_dir_name())
        self.assertEqual(summary.get_week_number(), '4')
        self.assertEqual(summary.get_year(), '2018')
        self.assertEqual(summary.get_csv(), '4,2018,0,0,0,0')

    def test_blastnfilter_summary_parse_summary_file_on_emptyfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            open(os.path.join(temp_dir,
                              BlastNFilterTask.SUMMARY_TXT), 'a').close()
            summary = BlastNFilterSummary(temp_dir)
            self.assertEqual(summary.get_csv(), '0,0,0,0,0,0')
        finally:
            shutil.rmtree(temp_dir)

    def test_blastnfilter_summary_parse_summary_file_on_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            os.mkdir(os.path.join(temp_dir,
                                  BlastNFilterTask.SUMMARY_TXT))
            summary = BlastNFilterSummary(temp_dir)
            self.assertEqual(summary.get_csv(), '0,0,0,0,0,0')
        finally:
            shutil.rmtree(temp_dir)

    def test_blastnfilter_summary_parse_summary_file_on_validfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            f = open(os.path.join(temp_dir,
                                  BlastNFilterTask.SUMMARY_TXT), 'w')
            f.write('INPUT SUMMARY\n')
            f.write('  entries:                             221\n')
            f.write('  complexes:                           178\n')
            f.write('  dockable complexes:                   95\n')
            f.write('  monomers:                            145\n')
            f.write('  dockable monomers:                    71\n')
            f.write('  multimers:                            76\n')
            f.write('  dockable multimers:                   24\n\n')

            f.write('FILTERING CRITERIA\n')
            f.write('  No. of query sequences           <=    1\n')
            f.write('  No. of dockable ligands           =    1\n')
            f.write('  Percent identity                 >=    0.95\n')
            f.write('  Percent Coverage                 >=    0.9\n')
            f.write('  No. of hit sequences             <=    4\n')
            f.write('  Structure determination method:        '
                    '  x-ray diffraction\n\n')

            f.write('OUTPUT SUMMARY\n')
            f.write('  Targets found:                        67\n')
            f.write('  Target: 5fz7|Sequences: 1|Hits: 94|Candidates: 17|'
                    'Elected:4|PDBids: 5fz7,5fyz,5a3p,5a1f\n')

            f.flush()
            f.close()

            summary = BlastNFilterSummary(temp_dir)
            self.assertEqual(summary.get_csv(), '0,0,178,95,71,67')
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
