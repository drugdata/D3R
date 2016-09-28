#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'


import unittest
import os
import tempfile
import shutil
from datetime import datetime
from dateutil.tz import tzutc
from dateutil.tz import tzlocal
from datetime import timedelta
from mock import Mock
from d3r.celpp.filetransfer import FtpFileTransfer


"""
test_dataimport
----------------------------------

Tests for `dataimport` module.
"""

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.dataimport import ImportRetryCountExceededError
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp import util

from tests.celpp import test_task


class TestDataImportTask(unittest.TestCase):
    def setUp(self):
        pass

    def get_total_seconds(self, td):
        return (td.microseconds +
                (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            self.assertEqual(task.get_uploadable_files(), [])

            task.create_dir()
            # test empty dir
            self.assertEqual(task.get_uploadable_files(), [])

            # test with only compinchi
            open(task.get_components_inchi_file(), 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(task.get_components_inchi_file())

            # test with crystal file
            open(task.get_crystalph_tsv(), 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(task.get_components_inchi_file())
            flist.index(task.get_crystalph_tsv())

            # test with nonpolymer file
            open(task.get_nonpolymer_tsv(), 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(task.get_components_inchi_file())
            flist.index(task.get_crystalph_tsv())
            flist.index(task.get_nonpolymer_tsv())

            # test with sequence file
            open(task.get_sequence_tsv(), 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 4)
            flist.index(task.get_components_inchi_file())
            flist.index(task.get_crystalph_tsv())
            flist.index(task.get_nonpolymer_tsv())
            flist.index(task.get_sequence_tsv())

        finally:
            shutil.rmtree(temp_dir)

    def test_update_status_from_filesystem(self):
        params = D3RParameters()
        task = DataImportTask(None, params)
        test_task.try_update_status_from_filesystem(self, task)

    def test_get_nonpolymer_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_nonpolymer_tsv(),
                         '/foo/' + task.get_dir_name() +
                         '/new_release_structure_nonpolymer.tsv')

    def test_get_sequence_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_sequence_tsv(),
                         '/foo/' + task.get_dir_name() +
                         '/new_release_structure_sequence_canonical.tsv')

    def test_get_crystalph_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_crystalph_tsv(),
                         '/foo/' + task.get_dir_name() +
                         '/new_release_crystallization_pH.tsv')

    def test_get_participant_list_csv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_participant_list_csv(),
                         '/foo/' + task.get_dir_name() +
                         '/' + DataImportTask.PARTICIPANT_LIST_CSV)

    def test_append_standard_to_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            task.append_standard_to_files()
            self.assertTrue(os.path.isfile(task.get_nonpolymer_tsv()))
            self.assertTrue(os.path.isfile(task.get_sequence_tsv()))
            self.assertTrue(os.path.isfile(task.get_crystalph_tsv()))

            # now do it again, but this time make the append fail
            # cause the nonpolymer_tsv is a directory
            os.unlink(task.get_nonpolymer_tsv())
            os.makedirs(task.get_nonpolymer_tsv())
            task.append_standard_to_files()

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run_with_complete_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            open(os.path.join(task.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(), None)
            self.assertEquals(task._can_run, False)

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run_does_not_exist_or_error(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)

            # no make blast db
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(),
                              'makeblastdb task has notfound status')
            self.assertEquals(task._can_run, False)

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()

            # make blast db failed
            err_file = os.path.join(make_blast.get_dir(),
                                    D3RTask.ERROR_FILE)
            open(err_file, 'a').close()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task.get_error(),
                              'makeblastdb task has error status')
            self.assertEquals(task._can_run, False)

            os.remove(err_file)

            # make blast db success
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            self.assertEquals(task.can_run(), True)
            self.assertEquals(task.get_error(), None)
            self.assertEquals(task._can_run, True)

            task.create_dir()
            open(os.path.join(task.get_dir(),
                              D3RTask.ERROR_FILE), 'a').close()
            self.assertEquals(task.can_run(), False)
            self.assertEquals(task._can_run, False)
            self.assertEquals(task.get_error(),
                              task.get_dir_name() + ' already exists and ' +
                              'status is ' + D3RTask.ERROR_STATUS)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_but_can_run_flag_is_false(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        task._can_run = False
        task.run()

    def test_run_pdbfileurl_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            task.run()
            self.assertEquals(task.get_error(), 'cannot download files' +
                              ' cause pdbfileurl not set')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_can_run_flag_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir
            task = DataImportTask(temp_dir, params)
            task._can_run = False
            task.run()
            self.assertEquals(task.get_error(), None)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_compinchi_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            task.run()
            self.assertEquals(task.get_error(), 'cannot download files' +
                              ' cause compinchi not set')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_nonpolymer_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            task.run()
            self.assertEquals(task.get_error(), 'Unable to download file ' +
                              'from ' + params.pdbfileurl + ' to ' +
                              task.get_nonpolymer_tsv())
        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_sequence_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            open(os.path.join(temp_dir,
                              task.NONPOLYMER_TSV), 'a').close()
            task.run()
            self.assertEquals(task.get_error(), 'Unable to download file ' +
                              'from ' + params.pdbfileurl + ' to ' +
                              task.get_sequence_tsv())

        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_crystalph_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            open(os.path.join(temp_dir,
                              task.NONPOLYMER_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.SEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.OLDSEQUENCE_TSV), 'a').close()
            task.run()
            self.assertEquals(task.get_error(), 'Unable to download file ' +
                              'from ' + params.pdbfileurl + ' to ' +
                              task.get_crystalph_tsv())

        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_compinchi_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            open(os.path.join(temp_dir,
                              task.NONPOLYMER_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.SEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.OLDSEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.CRYSTALPH_TSV), 'a').close()
            task.run()
            self.assertEquals(task.get_error(), 'Unable to download file ' +
                              'from ' + params.compinchi + ' to ' +
                              task.get_components_inchi_file())

        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_success_except_participant_download_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            task = DataImportTask(temp_dir, params)
            task._retrysleep = 0
            open(os.path.join(temp_dir,
                              task.NONPOLYMER_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.SEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.OLDSEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.CRYSTALPH_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.COMPINCHI_ICH), 'a').close()

            task.run()
            self.assertEquals(task.get_error(), None)

            # check line count is 1 now which indicates
            # standard was added
            self.assertEqual(util.get_file_line_count(
                task.get_nonpolymer_tsv()), 1)
            self.assertEqual(util.get_file_line_count(
                task.get_sequence_tsv()), 1)
            self.assertEqual(util.get_file_line_count(
                task.get_oldsequence_tsv()), 1)
            self.assertEqual(util.get_file_line_count(
                task.get_crystalph_tsv()), 1)

            self.assertTrue(task.get_email_log()
                            .startswith('\nWARNING: Unable to download'))
        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            fakeftp = FtpFileTransfer(None)
            mftp = D3RParameters()

            fakeftp.set_connection(mftp)
            fakeftp.set_remote_dir('/foo2')
            mftp.get = Mock()

            params = D3RParameters()
            params.pdbfileurl = 'file://' + temp_dir
            params.compinchi = 'file://' + temp_dir

            make_blast = MakeBlastDBTask(temp_dir, params)
            make_blast.create_dir()
            open(os.path.join(make_blast.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            task = DataImportTask(temp_dir, params)
            task.set_file_transfer(fakeftp)
            task._retrysleep = 0
            open(os.path.join(temp_dir,
                              task.NONPOLYMER_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.SEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.OLDSEQUENCE_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.CRYSTALPH_TSV), 'a').close()
            open(os.path.join(temp_dir,
                              task.COMPINCHI_ICH), 'a').close()

            task.run()
            self.assertEquals(task.get_error(), None)

            # check line count is 1 now which indicates
            # standard was added
            self.assertEqual(util.get_file_line_count(
                task.get_nonpolymer_tsv()), 1)
            self.assertEqual(util.get_file_line_count(
                task.get_sequence_tsv()), 1)
            self.assertEqual(util.get_file_line_count(
                task.get_oldsequence_tsv()), 1)
            self.assertEqual(util.get_file_line_count(
                task.get_crystalph_tsv()), 1)

            mftp.get.assert_called_with('/foo2/' +
                                        DataImportTask.PARTICIPANT_LIST_CSV,
                                        local=task.get_participant_list_csv())
        finally:
            shutil.rmtree(temp_dir)

    def test_wait_for_url_to_be_updated(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # skipimportwait is None
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task._wait_for_url_to_be_updated('foo')

            # skipimportwait is True
            params = D3RParameters()
            params.skipimportwait = True
            task = DataImportTask(temp_dir, params)
            task._wait_for_url_to_be_updated('foo')

            params = D3RParameters()
            params.skipimportwait = False
            params.importretry = 2
            params.importsleep = 0
            task = DataImportTask(temp_dir, params)

            fakefile = os.path.join(temp_dir, 'foo')

            f = open(fakefile, 'w')
            f.write('hi\n')
            f.flush()
            f.close()

            prev_friday = util.get_previous_friday_from_date(
                datetime.now(tzlocal()))

            thurs = prev_friday - timedelta(days=1)
            dse = thurs - datetime(1970, 01, 01, tzinfo=tzutc())
            secs_since_epoch = self.get_total_seconds(dse)
            os.utime(fakefile, (secs_since_epoch, secs_since_epoch))
            # test retry exceeded
            try:
                task._wait_for_url_to_be_updated('file://' + fakefile)
                self.fail('Expected ImportRetryCountExceededError')
            except ImportRetryCountExceededError:
                pass

            # test importsleep not defined
            params = D3RParameters()
            params.skipimportwait = False
            params.importretry = 2
            task.set_args(params)
            try:
                task._wait_for_url_to_be_updated('file://' + fakefile)
                self.fail('Expected ImportRetryCountExceededError')
            except ImportRetryCountExceededError:
                pass

            sat = prev_friday + timedelta(days=1)
            dse = sat - datetime(1970, 01, 01, tzinfo=tzutc())
            secs_since_epoch = self.get_total_seconds(dse)
            os.utime(fakefile, (secs_since_epoch, secs_since_epoch))

            task._wait_for_url_to_be_updated('file://' + fakefile)

        finally:
            shutil.rmtree(temp_dir)

    def test_download_files_with_wait_all_successful(self):
        temp_dir = tempfile.mkdtemp()
        try:
            prev_friday = util.get_previous_friday_from_date(
                datetime.now(tzlocal()))
            sat = prev_friday + timedelta(days=1)
            dse = sat - datetime(1970, 01, 01, tzinfo=tzutc())
            secs_since_epoch = self.get_total_seconds(dse)

            nonpoly = os.path.join(temp_dir,
                                   DataImportTask.NONPOLYMER_TSV)
            f = open(nonpoly, 'w')
            f.write('nonpoly\n')
            f.flush()
            f.close()
            os.utime(nonpoly, (secs_since_epoch, secs_since_epoch))

            seq = os.path.join(temp_dir,
                               DataImportTask.SEQUENCE_TSV)
            f = open(seq, 'w')
            f.write('seq\n')
            f.flush()
            f.close()
            os.utime(seq, (secs_since_epoch, secs_since_epoch))

            oldseq = os.path.join(temp_dir,
                                  DataImportTask.OLDSEQUENCE_TSV)
            f = open(oldseq, 'w')
            f.write('oldseq\n')
            f.flush()
            f.close()
            os.utime(oldseq, (secs_since_epoch, secs_since_epoch))

            crystal = os.path.join(temp_dir,
                                   DataImportTask.CRYSTALPH_TSV)
            f = open(crystal, 'w')
            f.write('crystal\n')
            f.flush()
            f.close()
            os.utime(crystal, (secs_since_epoch, secs_since_epoch))

            comp = os.path.join(temp_dir, DataImportTask.COMPINCHI_ICH)

            f = open(comp, 'w')
            f.write('comp\n')
            f.flush()
            f.close()

            params = D3RParameters()
            params.skipimportwait = False
            params.importretry = 2
            params.importsleep = 0
            params.compinchi = 'file://' + temp_dir
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            val = task._download_files('file://' + temp_dir)
            self.assertEqual(task.get_error(), None)
            self.assertTrue(val)

        finally:
            shutil.rmtree(temp_dir)

    def test_download_files_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            nonpoly = os.path.join(temp_dir,
                                   DataImportTask.NONPOLYMER_TSV)
            f = open(nonpoly, 'w')
            f.write('nonpoly\n')
            f.flush()
            f.close()

            seq = os.path.join(temp_dir,
                               DataImportTask.SEQUENCE_TSV)
            f = open(seq, 'w')
            f.write('seq\n')
            f.flush()
            f.close()

            oldseq = os.path.join(temp_dir,
                                  DataImportTask.OLDSEQUENCE_TSV)
            f = open(oldseq, 'w')
            f.write('oldseq\n')
            f.flush()
            f.close()

            crystal = os.path.join(temp_dir,
                                   DataImportTask.CRYSTALPH_TSV)
            f = open(crystal, 'w')
            f.write('crystal\n')
            f.flush()
            f.close()

            params = D3RParameters()
            params.skipimportwait = True
            params.compinchi = 'file://' + temp_dir
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            val = task._download_files('file://' + temp_dir)
            self.assertEqual(task.get_error(), 'Unable to download file from '
                                               'file://' + temp_dir + ' to ' +
                             task.get_components_inchi_file())
            self.assertFalse(val)

        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_from_crystalph_tsv_nonexistant_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            pdbid_set = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(len(pdbid_set), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_from_crystalph_tsv_empty_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            open(task.get_crystalph_tsv(), 'a').close()
            pdbid_set = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(os.path.isfile(task.get_crystalph_tsv()), True)
            self.assertEqual(len(pdbid_set), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_from_crystalph_tsv_no_header_invalid_lines(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('hello\nx\n')
            f.flush()
            f.close()
            pdbid_set = task.get_set_of_pdbid_from_crystalph_tsv()
            self.assertEqual(os.path.isfile(task.get_crystalph_tsv()), True)
            self.assertEqual(len(pdbid_set), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_from_crystalph_tsv_internal_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            open(task.get_crystalph_tsv(), 'a').close()
            os.chmod(task.get_crystalph_tsv(), 0o000)
            pdbid_set = task.get_set_of_pdbid_from_crystalph_tsv()
            self.assertEqual(len(pdbid_set), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_from_crystalph_tsv_valid_entries(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('PDB_ID  _exptl_crystal_grow.pH\n')
            f.write('4rfr    7.5\n')
            f.write('4X09\t6.5\n')
            f.write('4XET\t6.2\n')
            f.write('4XF1\t6.2\n')
            f.write('4XF3\t6.2\n')
            f.flush()
            f.close()
            pdbid_set = task.get_set_of_pdbid_from_crystalph_tsv()
            self.assertEqual(os.path.isfile(task.get_crystalph_tsv()), True)
            self.assertEqual(len(pdbid_set), 5)
            self.assertEqual('4RFR' in pdbid_set, True)
            self.assertEqual('4X09' in pdbid_set, True)
            self.assertEqual('4XET' in pdbid_set, True)
            self.assertEqual('4XF1' in pdbid_set, True)
            self.assertEqual('4XF3' in pdbid_set, True)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_from_crystalph_tsv_invalid_entries(self):
        """Header missing and 4rfr has 3 columns
        """
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('4X09\t6.5\n')
            f.write('4rfr 7.5 8\n')
            f.write('4XET\t6.2\n')
            f.write('4XF1\t6.2\n')
            f.write('4XF3\t6.2\n')
            f.flush()
            f.close()
            pdbid_set = task.get_set_of_pdbid_from_crystalph_tsv()
            self.assertEqual(os.path.isfile(task.get_crystalph_tsv()), True)
            self.assertEqual(len(pdbid_set), 3)
            self.assertEqual('4RFR' in pdbid_set, False)
            self.assertEqual('4X09' in pdbid_set, False)
            self.assertEqual('4XET' in pdbid_set, True)
            self.assertEqual('4XF1' in pdbid_set, True)
            self.assertEqual('4XF3' in pdbid_set, True)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres_no_seq(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('PDB_ID  _exptl_crystal_grow.pH\n')
            f.write('4X09\t6.5\n')
            f.write('4rfr\t8\n')
            f.write('4XET\t6.2\n')
            f.write('4XF1\t6.2\n')
            f.write('4XF3\t6.2\n')
            f.flush()
            f.close()
            pdbset = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres_no_tsv(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            makeblast = MakeBlastDBTask(temp_dir, params)
            makeblast.create_dir()
            f = open(makeblast.get_pdb_seqres_txt(), 'w')
            f.write('>101m_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVK'
                    'HLKTEAEMKASEDLKKHG\n')
            f.write('>102l_A mol:protein length:165  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKA'
                    'IGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAL'
                    'INMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRV'
                    'ITTFRTGTWDAYKNL\n')
            f.flush()
            f.close()

            pdbset = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres_empty_seq(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('PDB_ID  _exptl_crystal_grow.pH\n')
            f.write('4X09\t6.5\n')
            f.write('4rfr\t8\n')
            f.write('4XET\t6.2\n')
            f.write('4XF1\t6.2\n')
            f.write('4XF3\t6.2\n')
            f.flush()
            f.close()

            makeblast = MakeBlastDBTask(temp_dir, params)
            makeblast.create_dir()
            open(makeblast.get_pdb_seqres_txt(), 'a').close()

            pdbset = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres_no_hits(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('PDB_ID  _exptl_crystal_grow.pH\n')
            f.write('4X09\t6.5\n')
            f.write('4rfr\t8\n')
            f.write('4XET\t6.2\n')
            f.write('4XF1\t6.2\n')
            f.write('4XF3\t6.2\n')
            f.flush()
            f.close()

            makeblast = MakeBlastDBTask(temp_dir, params)
            makeblast.create_dir()
            f = open(makeblast.get_pdb_seqres_txt(), 'w')
            f.write('>101m_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVK'
                    'HLKTEAEMKASEDLKKHG\n')
            f.write('>102l_A mol:protein length:165  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKA'
                    'IGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAL'
                    'INMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRV'
                    'ITTFRTGTWDAYKNL\n')
            f.flush()
            f.close()

            pdbset = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres_w_hits(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_crystalph_tsv(), 'w')
            f.write('PDB_ID  _exptl_crystal_grow.pH\n')
            f.write('4X09\t6.5\n')
            f.write('4rfr\t8\n')
            f.write('4XET\t6.2\n')
            f.write('4XF1\t6.2\n')
            f.write('4XF3\t6.2\n')
            f.flush()
            f.close()

            makeblast = MakeBlastDBTask(temp_dir, params)
            makeblast.create_dir()
            f = open(makeblast.get_pdb_seqres_txt(), 'w')
            f.write('>4rfr_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVK'
                    'HLKTEAEMKASEDLKKHG\n')
            f.write('>102l_A mol:protein length:165  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKA'
                    'IGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAL'
                    'INMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRV'
                    'ITTFRTGTWDAYKNL\n')
            f.flush()
            f.close()

            pdbset = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(len(pdbset), 1)
            self.assertEqual('4RFR' in pdbset, True)
        finally:
            shutil.rmtree(temp_dir)

    def test_download_participant_list_no_filetransfer(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.create_dir()
            task._download_participant_list_csv()
            self.assertEqual(task.get_email_log(),
                             '\nWARNING: Unable to download '
                             'participant_list.csv which means '
                             'external users will NOT get evaluation '
                             'email : \'NoneType\' object has no attribute '
                             '\'connect\'\n')
        finally:
            shutil.rmtree(temp_dir)

    def test_download_participant_list_file_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileTransfer(None)
            mockftp = D3RParameters()
            mockftp.get = Mock()
            foo.set_remote_dir('/foo')
            foo.set_connection(mockftp)

            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.set_file_transfer(foo)
            task.create_dir()
            task._download_participant_list_csv()
            self.assertEqual(task.get_email_log(),
                             '\nWARNING: participant_list.csv not downloaded '
                             'which means external users will NOT get '
                             'evaluation email\n')
            mockftp.get\
                .assert_called_with('/foo/' +
                                    DataImportTask.PARTICIPANT_LIST_CSV,
                                    local=task.get_participant_list_csv())
        finally:
            shutil.rmtree(temp_dir)

    def test_download_participant_list_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileTransfer(None)
            mockftp = D3RParameters()
            mockftp.get = Mock()
            foo.set_remote_dir('/foo')
            foo.set_connection(mockftp)

            params = D3RParameters()
            task = DataImportTask(temp_dir, params)
            task.set_file_transfer(foo)
            task.create_dir()
            open(task.get_participant_list_csv(), 'a').close()
            task._download_participant_list_csv()
            self.assertEqual(task.get_email_log(), None)
            mockftp.get\
                .assert_called_with('/foo/' +
                                    DataImportTask.PARTICIPANT_LIST_CSV,
                                    local=task.get_participant_list_csv())
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
