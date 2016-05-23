#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import shutil
import os.path
from mock import Mock

from d3r.celpp.uploader import FtpFileUploader
"""
test_uploader
----------------------------------

Tests for `uploader` module.
"""


class MockFtp(object):
    """Dummy ftp connection
    """


class TestFtpFileUploader(unittest.TestCase):

    def setUp(self):
        pass

    def test_constructor(self):
        # test passing None which means we set everything manually
        foo = FtpFileUploader(None)

        self.assertEqual(foo.get_ftp_host(), None)
        self.assertEqual(foo.get_ftp_password(), None)
        self.assertEqual(foo.get_ftp_user(), None)
        self.assertEqual(foo.get_ftp_remote_dir(), '')
        self.assertEqual(foo.get_error_msg(), None)
        self.assertEqual(foo.get_ftp_remote_challenge_dir(), None)

        temp_dir = tempfile.mkdtemp()
        try:
            # test passing non existant config file
            try:
                foo = FtpFileUploader(os.path.join(temp_dir, 'doesnotexist'))
                self.fail('Expected IOError')
            except IOError:
                pass

            # test passing empty config file
            open(os.path.join(temp_dir, 'empty'), 'a').close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'empty'))
            self.assertEqual(foo.get_ftp_host(), None)
            self.assertEqual(foo.get_ftp_password(), None)
            self.assertEqual(foo.get_ftp_user(), None)
            self.assertEqual(foo.get_ftp_remote_dir(), '')
            self.assertEqual(foo.get_ftp_remote_challenge_dir(), None)

            # test getters and setters
            foo.set_connect_timeout(10)
            foo.set_ftp_connection('hi')
            foo.set_ftp_host('host')
            foo.set_ftp_password('pass')
            foo.set_ftp_remote_dir('/remote')
            foo.set_ftp_user('user')
            foo.set_ftp_remote_challenge_dir('/chall')

            self.assertEqual(foo.get_connect_timeout(), 10)
            self.assertEqual(foo.get_ftp_host(), 'host')
            self.assertEqual(foo.get_ftp_password(), 'pass')
            self.assertEqual(foo.get_ftp_remote_dir(), '/remote')
            self.assertEqual(foo.get_ftp_user(), 'user')
            self.assertEqual(foo.get_ftp_remote_challenge_dir(), '/chall')
        finally:
            shutil.rmtree(temp_dir)

    def test_parse_config(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test passing invalid config file
            f = open(os.path.join(temp_dir, 'invalid'), 'a')
            f.write('hello\nhow\nare')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'invalid'))
            self.assertEqual(foo.get_ftp_host(), None)
            self.assertEqual(foo.get_ftp_password(), None)
            self.assertEqual(foo.get_ftp_user(), None)
            self.assertEqual(foo.get_ftp_remote_dir(), '')
            self.assertEqual(foo.get_ftp_remote_challenge_dir(), None)

            # test passing partial config file
            f = open(os.path.join(temp_dir, 'partial'), 'a')
            f.write('#hello\nuser bob@bob.com\nare')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'partial'))
            self.assertEqual(foo.get_ftp_host(), None)
            self.assertEqual(foo.get_ftp_password(), None)
            self.assertEqual(foo.get_ftp_user(), 'bob@bob.com')
            self.assertEqual(foo.get_ftp_remote_dir(), '')
            self.assertEqual(foo.get_ftp_remote_challenge_dir(), None)

            # test passing valid config file
            f = open(os.path.join(temp_dir, 'valid'), 'a')
            f.write('host ftp.box.com\nuser bob@bob.com\npass 222\npath /foo')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_ftp_host(), 'ftp.box.com')
            self.assertEqual(foo.get_ftp_password(), '222')
            self.assertEqual(foo.get_ftp_user(), 'bob@bob.com')
            self.assertEqual(foo.get_ftp_remote_dir(), '/foo')
            self.assertEqual(foo.get_ftp_remote_challenge_dir(), None)

           # test passing valid config file with challengepath
            f = open(os.path.join(temp_dir, 'valid'), 'a')
            f.write('host ftp.box.com\nuser bob@bob.com\npass 222\npath /foo\n'
                    'challengepath /chall\n')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_ftp_host(), 'ftp.box.com')
            self.assertEqual(foo.get_ftp_password(), '222')
            self.assertEqual(foo.get_ftp_user(), 'bob@bob.com')
            self.assertEqual(foo.get_ftp_remote_dir(), '/foo')
            self.assertEqual(foo.get_ftp_remote_challenge_dir(), '/chall')

        finally:
            shutil.rmtree(temp_dir)

    def test_connect(self):
        # test where alt_ftp_con is None and we use ftpretty
        foo = FtpFileUploader(None)
        foo.set_ftp_host('doesnotexist')
        foo.set_ftp_user('user')
        foo.set_ftp_password('')
        foo.set_ftp_remote_dir('/remote')
        foo.set_connect_timeout(0)
        try:
            foo._connect()
            self.fail('expected exception')
        except:
            pass

        # test where alt_ftp_con is set
        foo = FtpFileUploader(None)
        foo.set_ftp_connection('hi')
        foo._connect()
        self.assertEqual(foo._ftp, 'hi')

    def test_disconnect(self):
        # test disconnect where _ftp is None
        foo = FtpFileUploader(None)
        foo._disconnect()

        # test disconnect where _ftp is set and
        # so is _alt_ftp_con
        foo = FtpFileUploader(None)
        foo.set_ftp_connection('hi')
        foo._connect()
        foo._disconnect()

        # test disconnect where _ftp is set and
        # _alt_ftp_con is None
        foo = FtpFileUploader(None)
        mockftp = MockFtp()
        mockftp.close = Mock(side_effect=Exception())
        foo.set_ftp_connection(mockftp)
        foo._connect()
        foo._alt_ftp_con = None
        foo._disconnect()
        mockftp.close.assert_any_call()

    def test_upload_file_direct_none_passed_as_file(self):
        foo = FtpFileUploader(None)

        self.assertEqual(foo.upload_file_direct(None,'/remote_dir',
                                                'remote_filename'), False)

        self.assertEqual(foo.get_error_msg(), 'File passed in is None')

    def test_upload_file_direct_file_is_not_a_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileUploader(None)
            noexist = os.path.join(temp_dir, 'noexist')
            self.assertEqual(foo.upload_file_direct(noexist,
                                   '/remote_dir', 'remote_filename'), False)

            self.assertEqual(foo.get_error_msg(), noexist +
                             ' is not a file')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_direct_remote_dir_is_none(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileUploader(None)
            afile = os.path.join(temp_dir, 'afile')
            open(afile, 'a').close()
            self.assertEqual(foo.upload_file_direct(afile, None,
                                                    'remote_filename'), False)
            self.assertEqual(foo.get_error_msg(), 'remote_dir is None')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_direct_remote_file_name_is_none(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileUploader(None)
            afile = os.path.join(temp_dir, 'afile')
            open(afile, 'a').close()
            self.assertEqual(foo.upload_file_direct(afile, '/foo', None),
                             False)
            self.assertEqual(foo.get_error_msg(), 'remote_file_name is None')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_direct_put_throws_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockftp = MockFtp()
            mockftp.put = Mock(side_effect=IOError('hi'))
            foo = FtpFileUploader(None)
            foo.set_ftp_connection(mockftp)
            afile = os.path.join(temp_dir, 'afile')
            open(afile, 'a').close()
            self.assertEqual(foo.upload_file_direct(afile, '/foo', 'name'),
                             False)
            self.assertEqual(foo.get_error_msg(), 'Unable to upload ' + afile +
                             ' to /foo/name : hi')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_direct_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockftp = MockFtp()
            mockftp.put = Mock(return_value=3)
            foo = FtpFileUploader(None)
            foo.set_ftp_connection(mockftp)
            afile = os.path.join(temp_dir, 'afile')
            open(afile, 'a').close()
            self.assertEqual(foo.upload_file_direct(afile, '/foo', 'name'),
                             True)
            self.assertEqual(foo.get_error_msg(), None)
            self.assertEqual(foo.get_upload_summary(), '1 (3 bytes) files '
                                                       'uploaded in 0 seconds'
                                                       ' to host Unset:')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_where_file_is_none(self):

        # test where file is None
        ftpspec = ["put"]
        mockftp = Mock(spec=ftpspec)
        foo = FtpFileUploader(None)
        foo.set_ftp_connection(mockftp)
        foo._connect()
        foo._upload_file(None)

    def test_upload_file_where_file_does_not_exist(self):
        temp_dir = tempfile.mkdtemp()
        try:
            ftpspec = ["put"]
            mockftp = Mock(spec=ftpspec)
            foo = FtpFileUploader(None)
            foo.set_ftp_connection(mockftp)
            foo._connect()
            foo._upload_file(os.path.join(temp_dir, 'nonexist'))
            self.assertEqual(mockftp.call_count, 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_that_raises_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:

            mockftp = MockFtp()
            mockftp.put = Mock(side_effect=IOError('hi'))
            foo = FtpFileUploader(None)
            foo.set_ftp_connection(mockftp)
            foo._connect()
            valid_file = os.path.join(temp_dir, 'file')
            f = open(valid_file, 'a')
            f.write('12')
            f.flush()
            f.close()
            try:
                foo._upload_file(valid_file)
                self.fail('Expected IOError')
            except IOError:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_on_valid_file(self):

        temp_dir = tempfile.mkdtemp()
        try:
            mockftp = MockFtp()
            mockftp.put = Mock(return_value=3)
            foo = FtpFileUploader(None)
            foo.set_ftp_connection(mockftp)
            foo.set_ftp_remote_dir('/remote')
            foo._connect()
            valid_file = os.path.join(temp_dir, 'file')
            f = open(valid_file, 'a')
            f.write('12')
            f.flush()
            f.close()
            foo._upload_file(valid_file)
            mockftp.put.assert_called_with(valid_file,
                                           os.path.normpath('/remote' +
                                                            valid_file))
            self.assertEqual(foo._bytes_transferred, 3)
            self.assertEqual(foo._files_transferred, 1)
            foo._upload_file(valid_file)
            self.assertEqual(foo._bytes_transferred, 6)
            self.assertEqual(foo._files_transferred, 2)
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_files_file_list_is_none(self):
        foo = FtpFileUploader(None)
        self.assertEqual(foo.upload_files(None), True)
        self.assertEqual(foo.get_error_msg(),
                         'List of files passed in was None')
        self.assertEqual(foo.get_upload_summary(),
                         'List of files passed in was None\n'
                         '0 (0 bytes) files uploaded in 0 '
                         'seconds to host Unset:')

    def test_upload_files_file_list_is_empty(self):
        foo = FtpFileUploader(None)
        self.assertEqual(foo.upload_files([]), True)
        self.assertEqual(foo.get_error_msg(),
                         'No files to upload')
        self.assertEqual(foo.get_upload_summary(),
                         'No files to upload\n'
                         '0 (0 bytes) files uploaded in 0 '
                         'seconds to host Unset:')

    def test_upload_files_connect_raises_exception(self):
        foo = FtpFileUploader(None)
        foo.set_ftp_host('ftp.box.com')
        foo.set_ftp_remote_dir('/hi')
        self.assertEqual(foo.upload_files(['hi']), False)
        self.assertEqual(foo.get_error_msg(), 'Unable to connect to ftp host')
        self.assertEqual(foo.get_upload_summary(),
                         'Unable to connect to ftp host\n'
                         '0 (0 bytes) files uploaded in 0 '
                         'seconds to host ftp.box.com:/hi')

    def test_upload_files_one_valid_file(self):
        mockftp = MockFtp()
        mockftp.put = Mock(return_value=3)
        mockftp.close = Mock(return_value=None)
        foo = FtpFileUploader(None)
        foo.set_ftp_connection(mockftp)

        foo.set_ftp_remote_dir('/remote')
        foo.set_ftp_host('hosty')
        temp_dir = tempfile.mkdtemp()
        try:
            valid_file = os.path.join(temp_dir, 'valid')
            f = open(valid_file, 'a')
            f.write('hi')
            f.flush()
            f.close()

            self.assertEqual(foo.upload_files([valid_file]), True)
            self.assertEqual(foo.get_error_msg(), None)
            self.assertEqual(foo.get_upload_summary(),
                             '1 (3 bytes) files uploaded in 0 '
                             'seconds to host hosty:/remote')
            mockftp.put.assert_called_with(valid_file,
                                           os.path.normpath('/remote' +
                                                            valid_file))
            mockftp.close.assert_not_called()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_files_multiple_valid_files(self):
        mockftp = MockFtp()
        mockftp.put = Mock(return_value=3)
        mockftp.close = Mock(return_value=None)
        foo = FtpFileUploader(None)
        foo.set_ftp_connection(mockftp)

        foo.set_ftp_remote_dir('/remote')
        foo.set_ftp_host('hosty')
        temp_dir = tempfile.mkdtemp()
        try:
            valid_file = os.path.join(temp_dir, 'valid')
            f = open(valid_file, 'a')
            f.write('hi')
            f.flush()
            f.close()
            afile = os.path.join(temp_dir, 'hi')
            open(afile, 'a').close()

            self.assertEqual(foo.upload_files([valid_file, afile]), True)
            self.assertEqual(foo.get_error_msg(), None)
            self.assertEqual(foo.get_upload_summary(),
                             '2 (6 bytes) files uploaded in 0 '
                             'seconds to host hosty:/remote')

            self.assertEqual(mockftp.put.call_count, 2)

            mockftp.close.assert_not_called()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_files_where_file_upload_raises_exception(self):
        mockftp = MockFtp()
        mockftp.put = Mock(side_effect=IOError('hi'))
        mockftp.close = Mock(return_value=None)
        foo = FtpFileUploader(None)
        foo.set_ftp_connection(mockftp)

        foo.set_ftp_remote_dir('/remote')
        foo.set_ftp_host('hosty')
        temp_dir = tempfile.mkdtemp()
        try:
            valid_file = os.path.join(temp_dir, 'valid')
            f = open(valid_file, 'a')
            f.write('hi')
            f.flush()
            f.close()

            self.assertEqual(foo.upload_files([valid_file]), False)
            self.assertEqual(foo.get_error_msg(), 'Error during upload')
            self.assertEqual(foo.get_upload_summary(),
                             'Error during upload\n'
                             '0 (0 bytes) files uploaded in 0 '
                             'seconds to host hosty:/remote')
            mockftp.put.assert_called_with(valid_file,
                                           os.path.normpath('/remote' +
                                                            valid_file))
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
