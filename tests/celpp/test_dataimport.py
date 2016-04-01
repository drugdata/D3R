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

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
