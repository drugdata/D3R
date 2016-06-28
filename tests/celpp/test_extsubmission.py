#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'

import unittest
import tempfile
import shutil
import os
from mock import Mock
from d3r.celpp.extsubmission import ExternalDataSubmissionFactory
from d3r.celpp.task import D3RParameters
from d3r.celpp import util

"""
test_extsubmission
--------------------------------

Tests for `extsubmission` module.
"""


class TestExternalSubmission(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_externaldatasubmissionfactory_noftpconfig(self):
        params = D3RParameters()
        fac = ExternalDataSubmissionFactory('/foo', params)
        self.assertEqual(fac.get_file_transfer(), None)

    def test_externaldatasubmissionfactory_ftpconfig_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            figfile = os.path.join(temp_dir,'foo')
            f = open(figfile, 'w')
            f.write('host blah.blah.com\n')
            f.write('user bob\n')
            f.write('pass ha\n')
            f.write('path /celppweekly/ha\n')
            f.write('challengepath /celppweekly/ha\n')
            f.flush()
            f.close()
            params = D3RParameters()
            params.ftpconfig = figfile
            fac = ExternalDataSubmissionFactory('/foo', params)
            self.assertTrue(fac.get_file_transfer() != None)
        finally:
            shutil.rmtree(temp_dir)

    def test_externaldatasubmissionfactory_getter_setter(self):
        params = D3RParameters()
        params.ha = 'hi'
        fac = ExternalDataSubmissionFactory('/foo', params)
        self.assertEqual(fac.get_file_transfer(), None)
        self.assertEqual(fac.get_path(), '/foo')
        self.assertEqual(fac.get_args().ha, 'hi')

        fac.set_file_transfer('yo')
        self.assertEqual(fac.get_file_transfer(), 'yo')

    def test_get_challenge_package_results_file_name(self):
        params = D3RParameters()

        # test under non existant dir
        fac = ExternalDataSubmissionFactory('/foo', params)
        self.assertEqual(fac._get_challenge_package_results_file_name('xxx'),
                         'celpp_week0_0_dockedresults_xxx.tar.gz')

        # test under 2016 week 40
        temp_dir = tempfile.mkdtemp()
        try:
            year = os.path.join(temp_dir, '2016')
            week = os.path.join(year, util.DATA_SET_WEEK_PREFIX + '40')
            os.makedirs(week)
            fac = ExternalDataSubmissionFactory(week, params)
            self.assertEqual(fac._get_challenge_package_results_file_name('xxx'),
                             'celpp_week40_2016_dockedresults_xxx.tar.gz')

        finally:
            shutil.rmtree(temp_dir)

    def test_get_submission_dirs_filetransfer_is_none(self):
        params = D3RParameters()
        try:
            fac = ExternalDataSubmissionFactory('/foo', params)
            fac._get_submission_dirs('ha')
            self.fail('Expected AttributeError')
        except AttributeError:
            pass

    def test_get_submission_dirs_list_is_none(self):
        params = D3RParameters()
        mockft = D3RParameters()
        mockft.list_dirs = Mock(return_value=None)
        fac = ExternalDataSubmissionFactory('/foo', params)
        fac.set_file_transfer(mockft)
        self.assertEqual(len(fac._get_submission_dirs('ha')), 0)
        mockft.list_dirs.assert_called_with('ha')

    def test_get_submissions_dirs_empty(self):
        params = D3RParameters()
        mockft = D3RParameters()
        mockft.list_dirs = Mock(return_value=[])
        fac = ExternalDataSubmissionFactory('/foo', params)
        fac.set_file_transfer(mockft)
        self.assertEqual(len(fac._get_submission_dirs('ha')), 0)
        mockft.list_dirs.assert_called_with('ha')

    def test_get_submission_dirs_one_dir(self):
        params = D3RParameters()
        mockft = D3RParameters()
        mockft.list_dirs = Mock(return_value=['somedir'])
        fac = ExternalDataSubmissionFactory('/foo', params)
        fac.set_file_transfer(mockft)
        dlist = fac._get_submission_dirs('ha')
        self.assertEqual(dlist[0], 'somedir')
        mockft.list_dirs.assert_called_with('ha')

    def test_get_submission_dirs_two_dirs(self):
        params = D3RParameters()
        mockft = D3RParameters()
        mockft.list_dirs = Mock(return_value=['somedir', 'blah'])
        fac = ExternalDataSubmissionFactory('/foo', params)
        fac.set_file_transfer(mockft)
        dlist = fac._get_submission_dirs('ha')
        self.assertEqual(dlist[0], 'somedir')
        self.assertEqual(dlist[1], 'blah')
        mockft.list_dirs.assert_called_with('ha')

    def test_get_challenge_data_package_file_raise_exception(self):
        params = D3RParameters()
        fac = ExternalDataSubmissionFactory('/foo', params)
        try:
            fac._get_challenge_data_package_file('/remote', 'dname')
            self.fail('Expected AttributeError')
        except AttributeError:
            pass

    def test_get_challenge_data_package_file_none_for_files(self):
        params = D3RParameters()
        mockft = D3RParameters()
        mockft.list_files = Mock(return_value=None)
        fac = ExternalDataSubmissionFactory('/foo', params)
        fac.set_file_transfer(mockft)
        dlist = fac._get_challenge_data_package_file('/remote', 'dname')
        self.assertEqual(dlist, None)
        mockft.list_files.assert_called_with('/remote/dname')

    def test_get_challenge_data_package_file_no_files(self):
        params = D3RParameters()
        mockft = D3RParameters()
        mockft.list_files = Mock(return_value=[])
        fac = ExternalDataSubmissionFactory('/foo', params)
        fac.set_file_transfer(mockft)
        dlist = fac._get_challenge_data_package_file('/remote', 'dname')
        self.assertEqual(dlist, None)
        mockft.list_files.assert_called_with('/remote/dname')

    def test_get_challenge_data_package_file_no_match(self):
        temp_dir = tempfile.mkdtemp()
        try:
            dr = ExternalDataSubmissionFactory.DOCKEDRESULTS
            year = os.path.join(temp_dir, '2016')
            week = os.path.join(year, 'dataset.week.13')
            os.makedirs(week)
            params = D3RParameters()
            mockft = D3RParameters()
            mockft.list_files = Mock(return_value=
                                     ['hi', 'celpp_week13_2015_'+
                                      dr +
                                      '_dname.tar.gz',
                                      'celpp_week13_2016' + dr +
                                      'yuck.tar.gz'])
            fac = ExternalDataSubmissionFactory(week, params)
            fac.set_file_transfer(mockft)
            dlist = fac._get_challenge_data_package_file('/remote', 'dname')
            self.assertEqual(dlist, None)
            mockft.list_files.assert_called_with('/remote/dname')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_challenge_data_package_file_found_match(self):
        temp_dir = tempfile.mkdtemp()
        try:
            dr = ExternalDataSubmissionFactory.DOCKEDRESULTS
            year = os.path.join(temp_dir, '2016')
            week = os.path.join(year, 'dataset.week.13')
            os.makedirs(week)
            params = D3RParameters()
            mockft = D3RParameters()
            mockft.list_files = Mock(return_value=
                                     ['hi', 'celpp_week13_2015' + dr +
                                      'dname.tar.gz',
                                      'celpp_week13_2016' + dr +
                                      'yuck.tar.gz'])
            fac = ExternalDataSubmissionFactory(week, params)
            fac.set_file_transfer(mockft)
            dlist = fac._get_challenge_data_package_file('/remote', 'yuck')
            self.assertEqual(dlist,
                             '/remote/yuck/celpp_week13_2016' + dr +
                             'yuck.tar.gz')
            mockft.list_files.assert_called_with('/remote/yuck')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_external_data_submissions_raise_exception(self):
        params = D3RParameters()
        fac = ExternalDataSubmissionFactory('/foo', params)
        tlist = fac.get_external_data_submissions()
        self.assertEqual(len(tlist), 0)

    def test_get_external_data_submissions_none_subdir(self):
        params = D3RParameters()
        fac = ExternalDataSubmissionFactory('/foo', params)
        mockft = D3RParameters()
        mockft.list_dirs = Mock(return_value=None)
        fac.set_file_transfer(mockft)
        tlist = fac.get_external_data_submissions()
        self.assertEqual(len(tlist), 0)

    def test_get_external_data_submissions_no_dirs(self):
        params = D3RParameters()
        fac = ExternalDataSubmissionFactory('/foo', params)
        mockft = D3RParameters()
        mockft.list_dirs = Mock(return_value=[])
        fac.set_file_transfer(mockft)
        tlist = fac.get_external_data_submissions()
        self.assertEqual(len(tlist), 0)

    def test_get_external_data_submissions_no_chall_packages(self):
        temp_dir = tempfile.mkdtemp()
        try:
            year = os.path.join(temp_dir, '2016')
            week = os.path.join(year, 'dataset.week.13')
            os.makedirs(week)
            dr = ExternalDataSubmissionFactory.DOCKEDRESULTS
            params = D3RParameters()
            fac = ExternalDataSubmissionFactory(week, params)
            mockft = D3RParameters()
            mockft.connect = Mock()
            mockft.disconnect = Mock()
            mockft.get_ftp_remote_submission_dir = Mock(return_value=
                                                        '/remote')
            mockft.list_dirs = Mock(return_value=['yo'])
            mockft.list_files = Mock(return_value=['hi',
                                                    'celpp_week13_2015' + dr +
                                                    'dname.tar.gz',
                                                    'celpp_week13_2016' + dr +
                                                    'yuck.tar.gz'])
            fac.set_file_transfer(mockft)
            tlist = fac.get_external_data_submissions()
            self.assertEqual(len(tlist), 0)
            mockft.list_dirs.assert_called_with('/remote')
            mockft.list_files.assert_called_with('/remote/yo')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_external_data_submissions_one_extdatasubmission(self):
        temp_dir = tempfile.mkdtemp()
        try:
            year = os.path.join(temp_dir, '2017')
            week = os.path.join(year, 'dataset.week.13')
            os.makedirs(week)
            dr = ExternalDataSubmissionFactory.DOCKEDRESULTS
            params = D3RParameters()
            fac = ExternalDataSubmissionFactory(week, params)
            mockft = D3RParameters()
            mockft.connect = Mock()
            mockft.disconnect = Mock()
            mockft.get_ftp_remote_challenge_dir = Mock(return_value='/chall')
            mockft.delete_file = Mock(return_value=False)

            mockft.get_ftp_remote_submission_dir = Mock(return_value=
                                                        '/remote')
            mockft.list_dirs = Mock(return_value=['yo'])
            mockft.list_files = Mock(return_value=['hi',
                                                    'celpp_week13_2017' + dr +
                                                    'yo.tar.gz',
                                                    'celpp_week13_2016' + dr +
                                                    'yuck.tar.gz'])
            fac.set_file_transfer(mockft)
            tlist = fac.get_external_data_submissions()
            self.assertEqual(len(tlist), 1)
            self.assertEqual(tlist[0].get_name(), 'yo.extsubmission')
            self.assertEqual(tlist[0].get_remote_challenge_data_package(),
                             '/remote/yo/celpp_week13_2017' + dr + 'yo.tar.gz')

            mockft.list_dirs.assert_called_with('/remote')
            mockft.list_files.assert_called_with('/remote/yo')
            mockft.delete_file.assert_called_with('/chall/latest.txt')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_external_data_submissions_two_extdatasubmission(self):
        temp_dir = tempfile.mkdtemp()
        try:
            year = os.path.join(temp_dir, '2017')
            week = os.path.join(year, 'dataset.week.13')
            os.makedirs(week)
            dr = ExternalDataSubmissionFactory.DOCKEDRESULTS
            params = D3RParameters()
            fac = ExternalDataSubmissionFactory(week, params)
            mockft = D3RParameters()
            mockft.connect = Mock()
            mockft.delete_file = Mock(return_value=True)
            mockft.get_ftp_remote_challenge_dir = Mock(return_value='/chall')
            mockft.get_ftp_remote_submission_dir = Mock(return_value=
                                                        '/remote')
            mockft.list_dirs = Mock(return_value=['yo', 'yuck'])
            mockft.list_files = Mock(return_value=['hi',
                                                    'celpp_week13_2017' + dr +
                                                    'yo.tar.gz',
                                                    'celpp_week13_2017' + dr +
                                                    'yuck.tar.gz'])
            fac.set_file_transfer(mockft)
            tlist = fac.get_external_data_submissions()
            self.assertEqual(len(tlist), 2)
            self.assertEqual(tlist[0].get_name(), 'yo.extsubmission')
            self.assertEqual(tlist[0].get_remote_challenge_data_package(),
                             '/remote/yo/celpp_week13_2017' + dr + 'yo.tar.gz')
            self.assertEqual(tlist[1].get_name(), 'yuck.extsubmission')
            self.assertEqual(tlist[1].get_remote_challenge_data_package(),
                             '/remote/yuck/celpp_week13_2017' + dr + 'yuck.tar.gz')
            mockft.list_dirs.assert_called_with('/remote')
            mockft.list_files.assert_called_with('/remote/yuck')
            mockft.delete_file.assert_called_with('/chall/latest.txt')
        finally:
            shutil.rmtree(temp_dir)


if __name__ == '__main__':
    unittest.main()
