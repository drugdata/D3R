#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'


import unittest
import os
import tempfile
import shutil

"""
test_dataimport
----------------------------------

Tests for `dataimport` module.
"""

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp import util
from tests.celpp import test_task


class TestDataImportTask(unittest.TestCase):
    def setUp(self):
        pass

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
                         '/new_release_structure_sequence.tsv')

    def test_get_crystalph_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_crystalph_tsv(),
                         '/foo/' + task.get_dir_name() +
                         '/new_release_crystallization_pH.tsv')

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
                              task.CRYSTALPH_TSV), 'a').close()
            task.run()
            self.assertEquals(task.get_error(), 'Unable to download file ' +
                              'from ' + params.compinchi + ' to ' +
                              task.get_components_inchi_file())

        finally:
            shutil.rmtree(temp_dir)

    def test_run_all_success(self):
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
                task.get_crystalph_tsv()), 1)
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
        except:
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
        except:
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
            pdbid_set = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(os.path.isfile(task.get_crystalph_tsv()), True)
            self.assertEqual(len(pdbid_set), 5)
            self.assertEqual('4RFR' in pdbid_set, True)
            self.assertEqual('4X09' in pdbid_set, True)
            self.assertEqual('4XET' in pdbid_set, True)
            self.assertEqual('4XF1' in pdbid_set, True)
            self.assertEqual('4XF3' in pdbid_set, True)
        except:
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
            pdbid_set = task.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            self.assertEqual(os.path.isfile(task.get_crystalph_tsv()), True)
            self.assertEqual(len(pdbid_set), 3)
            self.assertEqual('4RFR' in pdbid_set, False)
            self.assertEqual('4X09' in pdbid_set, False)
            self.assertEqual('4XET' in pdbid_set, True)
            self.assertEqual('4XF1' in pdbid_set, True)
            self.assertEqual('4XF3' in pdbid_set, True)
        except:
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


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
