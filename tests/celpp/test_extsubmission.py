#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'

import unittest
import tempfile
import shutil
import os
import tarfile
from mock import Mock

from d3r.celpp.extsubmission import ExternalDataSubmissionFactory
from d3r.celpp.extsubmission import ExternalDataSubmissionTask
from d3r.celpp.extsubmission import ChallengePackageDownloadError
from d3r.celpp.task import D3RParameters
from d3r.celpp import util
from d3r.celpp.task import D3RTask

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
            figfile = os.path.join(temp_dir, 'foo')
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
            self.assertTrue(fac.get_file_transfer() is not None)
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
            rfname = fac._get_challenge_package_results_file_name('xxx')
            self.assertEqual(rfname,
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
            mockft.list_files = Mock(return_value=['hi', 'celpp_week13_2015_' +
                                                   dr + '_dname.tar.gz',
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
            mockft.list_files = Mock(return_value=['hi', 'celpp_week13_2015' +
                                                   dr + 'dname.tar.gz',
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
            mockft.get_remote_submission_dir = Mock(return_value='/remote')
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
            mockft.get_remote_challenge_dir = Mock(return_value='/chall')
            mockft.delete_file = Mock(return_value=False)

            mockft.get_remote_submission_dir = Mock(return_value='/remote')
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
            mockft.get_remote_challenge_dir = Mock(return_value='/chall')
            mockft.get_remote_submission_dir = Mock(return_value='/remote')
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
                             '/remote/yuck/celpp_week13_2017' + dr +
                             'yuck.tar.gz')
            mockft.list_dirs.assert_called_with('/remote')
            mockft.list_files.assert_called_with('/remote/yuck')
            mockft.delete_file.assert_called_with('/chall/latest.txt')
        finally:
            shutil.rmtree(temp_dir)

    def test_externaltask_get_set_remote_challenge_data_package(self):
        params = D3RParameters()
        task = ExternalDataSubmissionTask('/foo', 'name', '/remote', params)
        self.assertEqual(task.get_remote_challenge_data_package(), '/remote')
        task.set_remote_challenge_data_package('/blah')
        self.assertEqual(task.get_remote_challenge_data_package(), '/blah')

    def test_externaltask_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ExternalDataSubmissionTask(temp_dir, 'name', '/remote',
                                              params)
            self.assertTrue(task.can_run())

            task.create_dir()
            self.assertFalse(task.can_run())
            self.assertEqual(task.get_error(), 'stage.6.name.extsubmission '
                                               'already exists and status is '
                                               'unknown')
            open(os.path.join(task.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            self.assertFalse(task.can_run())
            self.assertEqual(task.get_error(), None)
        finally:
            shutil.rmtree(temp_dir)

    def test_is_tarmember_safe(self):
        params = D3RParameters()
        task = ExternalDataSubmissionTask('/foo', 'name',
                                          '/remote/celpp_week12_2016_'
                                          'dockedresults_name.tar.gz', params)
        # path starts with /
        val = task._is_tarmembername_safe('/celpp_week12', 'celpp_week12')
        self.assertFalse(val)
        self.assertEqual(task.get_email_log(), 'Skipping, path starts with '
                                               '/ : /celpp_week12\n')
        task.set_email_log('')
        # path does not start with chall_name
        val = task._is_tarmembername_safe('celpp', 'celpp_week12')
        self.assertFalse(val)
        self.assertEqual(task.get_email_log(), 'Skipping, path does not'
                                               ' conform: celpp\n')
        task.set_email_log('')

        # path has .. in it
        val = task._is_tarmembername_safe('celpp/../yo', 'celpp')
        self.assertFalse(val)
        self.assertEqual(task.get_email_log(), 'Skipping, found .. in path:'
                                               ' celpp/../yo\n')
        task.set_email_log('')
        chall_name = 'celpp_week12_2016_dockedresults_name'
        # valid path
        val = task._is_tarmembername_safe(chall_name + '/5xjv/yo.pdb',
                                          chall_name)
        self.assertTrue(val)
        self.assertEqual(task.get_email_log(), '')

    def test_untar_challenge_data_package_no_tarfile(self):

        params = D3RParameters()
        task = ExternalDataSubmissionTask('/foo', 'name',
                                          '/remote/hello.tar.gz', params)
        try:
            task._untar_challenge_data_package('hello')
            self.fail('Expected IOError')
        except IOError:
            pass

    def create_good_tarfile(self, temp_dir, nameprefix):
        tfiledir = os.path.join(temp_dir, nameprefix)
        tfile = tfiledir + '.tar.gz'

        os.makedirs(tfiledir)
        candidateone = '5juw'
        candidateonedir = os.path.join(tfiledir, candidateone)
        os.makedirs(candidateonedir)

        mypdb = 'apo-5juw_2eb2_docked.pdb'
        mypdbfile = os.path.join(candidateonedir, mypdb)
        f = open(mypdbfile, 'w')
        f.write('mypdb')
        f.flush()
        f.close()

        mymol = 'apo-5juw_2eb2_docked.mol'
        mymolfile = os.path.join(candidateonedir, mymol)
        f = open(mymolfile, 'w')
        f.write('mymol')
        f.flush()
        f.close()

        tar = tarfile.open(tfile, 'w:gz')
        tar.add(candidateonedir, arcname=nameprefix + '/' + candidateone)
        tar.close()
        shutil.rmtree(tfiledir)

        return tfile

    def create_bad_tarfile_dotdot_path(self, temp_dir, nameprefix):
        tfiledir = os.path.join(temp_dir, nameprefix)
        tfile = tfiledir + '.tar.gz'

        os.makedirs(tfiledir)
        candidateone = '5juw'
        candidateonedir = os.path.join(tfiledir, candidateone)
        os.makedirs(candidateonedir)

        mypdb = 'apo-5juw_2eb2_docked.pdb'
        mypdbfile = os.path.join(candidateonedir, mypdb)
        f = open(mypdbfile, 'w')
        f.write('mypdb')
        f.flush()
        f.close()

        mymol = 'apo-5juw_2eb2_..docked.mol'
        mymolfile = os.path.join(candidateonedir, mymol)
        f = open(mymolfile, 'w')
        f.write('mymol')
        f.flush()
        f.close()

        # add in a symlink for fun
        os.symlink(mymol, os.path.join(candidateonedir, 'hello'))

        tar = tarfile.open(tfile, 'w:gz')
        tar.add(candidateonedir, arcname='/' + nameprefix + '/' + candidateone)
        tar.close()
        shutil.rmtree(tfiledir)

        return tfile

    def test_untar_challenge_data_package_good_tarfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ExternalDataSubmissionTask(temp_dir, 'hello',
                                              'hi', params)
            task.create_dir()
            self.create_good_tarfile(task.get_dir(), 'hello')
            c = task._untar_challenge_data_package('hello.tar.gz')
            self.assertEqual(c, 'hello')
            self.assertEqual(task.get_email_log(), None)
            candidate = os.path.join(task.get_dir(), 'hello', '5juw')
            self.assertTrue(os.path.isdir(candidate))
            mypdb = os.path.join(candidate,
                                 'apo-5juw_2eb2_docked.pdb')
            mymol = os.path.join(candidate,
                                 'apo-5juw_2eb2_docked.mol')
            self.assertTrue(os.path.isfile(mypdb))
            self.assertTrue(os.path.isfile(mymol))
        finally:
            shutil.rmtree(temp_dir)

    def test_untar_challenge_data_package_dotdot_tarfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ExternalDataSubmissionTask(temp_dir, 'hello',
                                              'hi', params)
            task.create_dir()
            self.create_bad_tarfile_dotdot_path(task.get_dir(),
                                                'hello')
            c = task._untar_challenge_data_package('hello.tar.gz')
            self.assertEqual(c, 'hello')
            self.assertTrue('Skipping, found .. in '
                            'path: hello/5juw/apo-'
                            '5juw_2eb2_..docked.mol\n' in task.get_email_log())

            self.assertTrue('Ignoring non dir/file '
                            'entry in tar 5juw/'
                            'hello\n' in task.get_email_log())

            candidate = os.path.join(task.get_dir(), 'hello', '5juw')
            self.assertTrue(os.path.isdir(candidate))
            mypdb = os.path.join(candidate,
                                 'apo-5juw_2eb2_docked.pdb')
            mymol = os.path.join(candidate,
                                 'apo-5juw_2eb2.._docked.mol')
            self.assertTrue(os.path.isfile(mypdb))
            self.assertFalse(os.path.isfile(mymol))
        finally:
            shutil.rmtree(temp_dir)

    def test_untar_challenge_data_package_wrong_prefix(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              'hi', params)
            task.create_dir()
            tfile = self.create_good_tarfile(task.get_dir(),
                                             'yo')
            shutil.move(tfile, os.path.join(task.get_dir(),
                                            'byte.tar.gz'))
            c = task._untar_challenge_data_package('byte.tar.gz')
            self.assertEqual(c, 'byte')

            candidate = os.path.join(task.get_dir(), 'hello', '5juw')
            self.assertFalse(os.path.isdir(candidate))
            mypdb = os.path.join(candidate,
                                 'apo-5juw_2eb2_docked.pdb')
            mymol = os.path.join(candidate,
                                 'apo-5juw_2eb2.._docked.mol')
            self.assertFalse(os.path.isfile(mypdb))
            self.assertFalse(os.path.isfile(mymol))
        finally:
            shutil.rmtree(temp_dir)

    def test_move_challenge_data_package_into_task_dir_raises_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              'hi', params)
            task.create_dir()
            try:
                task._move_challenge_data_package_into_task_dir('foo')
                self.fail('expected IOError')
            except OSError:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_move_challenge_data_package_into_task_dir_valid(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              'hi', params)
            task.create_dir()
            basedir = os.path.join(task.get_dir(), 'foo')
            os.makedirs(basedir)
            can1 = os.path.join(basedir, '5juw')
            os.makedirs(can1)
            pdb1 = os.path.join(can1, 'apo.pdb')
            open(pdb1, 'a').close()
            mol1 = os.path.join(can1, 'apo.mol')
            open(mol1, 'a').close()

            can2 = os.path.join(basedir, 'abcd')
            os.makedirs(can2)
            pdb2 = os.path.join(can2, 'apo.pdb')
            open(pdb2, 'a').close()
            mol2 = os.path.join(can2, 'apo.mol')
            open(mol2, 'a').close()

            task._move_challenge_data_package_into_task_dir('foo')
            can3 = os.path.join(task.get_dir(), '5juw')
            can4 = os.path.join(task.get_dir(), 'abcd')
            self.assertTrue(os.path.isdir(can3))
            self.assertTrue(os.path.isdir(can4))
            self.assertTrue(os.path.isfile(os.path.join(can3, 'apo.pdb')))
            self.assertTrue(os.path.isfile(os.path.join(can3, 'apo.mol')))
            self.assertTrue(os.path.isfile(os.path.join(can4, 'apo.pdb')))
            self.assertTrue(os.path.isfile(os.path.join(can4, 'apo.mol')))

        finally:
            shutil.rmtree(temp_dir)

    def test_get_summary_of_docked_results(self):
        params = D3RParameters()
        task = ExternalDataSubmissionTask('/foo', 'yo',
                                          'hi', params)
        self.assertEqual(task._get_summary_of_docked_results(), '')

    def test_download_remote_challenge_data_package_with_retry_retry_neg(self):
        params = D3RParameters()
        task = ExternalDataSubmissionTask('/foo', 'yo',
                                          'hi', params)
        task.set_download_max_retry_count(-1)
        task.set_download_retry_sleep(0)
        try:
            task._download_remote_challenge_data_package_with_retry()
            self.fail('Expected ChallengePackageDownloadError')
        except ChallengePackageDownloadError as e:
            self.assertEqual(str(e), 'Unable to download hi')

    def test_download_remote_challenge_data_package_with_retry_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp1_2017_dockedresults_yo.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            task.set_download_retry_sleep(0)
            task.set_download_max_retry_count(1)
            mockft = D3RParameters()
            mockft.connect = Mock(side_effect=IOError('error'))
            mockft.download_file = Mock(return_value=True)
            mockft.disconnect = Mock(return_value=None)
            task.set_file_transfer(mockft)
            try:
                task._download_remote_challenge_data_package_with_retry()
            except ChallengePackageDownloadError as e:
                self.assertEqual(str(e), 'Unable to download ' +
                                 'celpp1_2017_dockedresults_yo.tar.gz')
        finally:
            shutil.rmtree(temp_dir)

    def test_download_remote_challenge_data_package_with_retry_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp1_2017_dockedresults_yo.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            task.set_download_retry_sleep(0)
            task.set_download_max_retry_count(1)
            mockft = D3RParameters()
            mockft.connect = Mock()
            mockft.download_file = Mock(return_value=False)
            mockft.disconnect = Mock(return_value=None)
            task.set_file_transfer(mockft)
            try:
                task._download_remote_challenge_data_package_with_retry()
            except ChallengePackageDownloadError as e:
                self.assertEqual(str(e), 'Unable to download ' +
                                 'celpp1_2017_dockedresults_yo.tar.gz')
        finally:
            shutil.rmtree(temp_dir)

    def test_download_remote_challenge_data_package_with_1st_retry_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp1_2017_dockedresults_yo.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            task.set_download_retry_sleep(0)
            task.set_download_max_retry_count(3)
            mockft = D3RParameters()
            mockft.connect = Mock()
            mockft.download_file = Mock()
            mockft.download_file.side_effect = [False, True]
            mockft.disconnect = Mock(return_value=None)
            task.set_file_transfer(mockft)
            c = task._download_remote_challenge_data_package_with_retry()
            self.assertEqual(c, 'celpp1_2017_dockedresults_yo.tar.gz')
        finally:
            shutil.rmtree(temp_dir)

    def test__download_remote_challenge_data_package_raise_exception(self):
        params = D3RParameters()
        task = ExternalDataSubmissionTask('/foo', 'yo',
                                          'hi', params)
        try:
            task._download_remote_challenge_data_package('foo')
            self.fail('Expected attribute error cause file transfer'
                      'object is not set')
        except AttributeError:

            pass

    def test__download_remote_challenge_data_package_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp2_3.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            task.create_dir()
            mockft = D3RParameters()
            mockft.connect = Mock(return_value=None)
            mockft.download_file = Mock(return_value=True)
            mockft.disconnect = Mock(return_value=None)
            task.set_file_transfer(mockft)
            val = task._download_remote_challenge_data_package('celpp2_3.'
                                                               'tar.gz')
            self.assertTrue(val)
            localfile = os.path.join(task.get_dir(), 'celpp2_3.tar.gz')
            mockft.download_file.assert_called_with(pkg, localfile)
            mockft.connect.assert_called_with()
            mockft.disconnect.assert_called_with()
        finally:
            shutil.rmtree(temp_dir)

    def test_run_unable_to_run_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp2_3.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            task.create_dir()
            task.run()
            self.assertEqual(task.get_email_log(), None)
            self.assertEqual(task.get_error(), 'stage.6.yo.extsubmission '
                                               'already exists and status '
                                               'is unknown')
        finally:
            shutil.rmtree(temp_dir)

    def test_runtask_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp1_2017_dockedresults_yo.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            mockft = D3RParameters()
            mockft.connect = Mock(return_value=None)
            mockft.download_file = Mock(return_value=True)
            mockft.disconnect = Mock(return_value=None)
            task.set_file_transfer(mockft)
            task.create_dir()
            self.create_good_tarfile(task.get_dir(),
                                     'celpp1_2017_dockedresults_yo')
            task._runtask()
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_email_log(), '')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fail(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            pkg = '/remote/celpp1_2017_dockedresults_yo.tar.gz'
            task = ExternalDataSubmissionTask(temp_dir, 'yo',
                                              pkg, params)
            mockft = D3RParameters()
            mockft.connect = Mock(side_effect=IOError('error'))
            mockft.download_file = Mock(return_value=True)
            mockft.disconnect = Mock(return_value=None)
            task.set_file_transfer(mockft)
            task.run()
            self.assertEqual(task.get_error(), 'Caught exception Unable to '
                                               'download celpp1_2017_docked'
                                               'results_yo.tar.gz')
            self.assertEqual(task.get_email_log(), 'Try # 1 of 3 to download '
                                                   'celpp1_2017_dockedresults'
                                                   '_yo.tar.gz\nTry # 2 of 3 '
                                                   'to download celpp1_2017_d'
                                                   'ockedresults_yo.tar.gz\n'
                                                   'Try # 3 of 3 to download '
                                                   'celpp1_2017_dockedresults'
                                                   '_yo.tar.gz\n')

        finally:
            shutil.rmtree(temp_dir)


if __name__ == '__main__':
    unittest.main()
