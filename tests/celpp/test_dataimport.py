#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'


import unittest

"""
test_dataimport
----------------------------------

Tests for `dataimport` module.
"""

from d3r.celpp.task import D3RParameters
from d3r.celpp.dataimport import DataImportTask

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
                         '/foo/stage.1.dataimport' +
                         '/new_release_structure_nonpolymer.tsv')

    def test_get_sequence_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_sequence_tsv(),
                         '/foo/stage.1.dataimport/' +
                         'new_release_structure_sequence.tsv')

    def test_get_crystalph_tsv(self):
        params = D3RParameters()
        task = DataImportTask('/foo', params)
        self.assertEqual(task.get_crystalph_tsv(),
                         '/foo/stage.1.dataimport/' +
                         'new_release_crystallization_pH.tsv')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
