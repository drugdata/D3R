#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_compinchidownload
--------------------------------

Tests for `compinchidownload` module.
"""

import shutil
import os

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.compinchidownload import CompInchiDownloadTask


class TestCompInchiDownloadTask(unittest.TestCase):

    def setUp(self):
        pass

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
