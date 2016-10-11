#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import shutil
import os.path
from mock import Mock

from d3r.celpp.filetransfer import FtpFileTransfer

"""
test_uploader
----------------------------------

Tests for `uploader` module.
"""


class MockFtp(object):
    """Dummy ftp connection
    """


class TestFtpFileTransfer(unittest.TestCase):

    def setUp(self):
        pass

    def test_constructor(self):
        # test passing None which means we set everything manually
        foo = FtpFileTransfer(None)

        self.assertEqual(foo.get_host(), None)
        self.assertEqual(foo.get_password(), None)
        self.assertEqual(foo.get_user(), None)
        self.assertEqual(foo.get_remote_dir(), '')
        self.assertEqual(foo.get_error_msg(), None)
        self.assertEqual(foo.get_remote_challenge_dir(), None)

        temp_dir = tempfile.mkdtemp()
        try:
            # test passing non existant config file
            try:
                foo = FtpFileTransfer(os.path.join(temp_dir, 'doesnotexist'))
                self.fail('Expected IOError')
            except IOError:
                pass

            # test passing empty config file
            open(os.path.join(temp_dir, 'empty'), 'a').close()
            foo = FtpFileTransfer(os.path.join(temp_dir, 'empty'))
            self.assertEqual(foo.get_host(), None)
            self.assertEqual(foo.get_password(), None)
            self.assertEqual(foo.get_user(), None)
            self.assertEqual(foo.get_remote_dir(), '')
            self.assertEqual(foo.get_remote_challenge_dir(), None)
            self.assertEqual(foo.get_contestant_id(), None)

            # test getters and setters
            foo.set_connect_timeout(10)
            foo.set_connection('hi')
            foo.set_host('host')
            foo.set_password('pass')
            foo.set_remote_dir('/remote')
            foo.set_user('user')
            foo.set_remote_challenge_dir('/chall')
            foo.set_contestant_id(12345)

            self.assertEqual(foo.get_connect_timeout(), 10)
            self.assertEqual(foo.get_host(), 'host')
            self.assertEqual(foo.get_password(), 'pass')
            self.assertEqual(foo.get_remote_dir(), '/remote')
            self.assertEqual(foo.get_user(), 'user')
            self.assertEqual(foo.get_remote_challenge_dir(), '/chall')
            self.assertEqual(foo.get_contestant_id(), '12345')

            foo.set_contestant_id(None)
            self.assertEqual(foo.get_contestant_id(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_connect(self):
        # test where alt_ftp_con is None and we use ftpretty
        foo = FtpFileTransfer(None)
        foo.set_host('doesnotexist')
        foo.set_user('user')
        foo.set_password('')
        foo.set_remote_dir('/remote')
        foo.set_connect_timeout(0)
        self.assertFalse(foo.connect())

        # test where alt_ftp_con is set
        foo = FtpFileTransfer(None)
        foo.set_connection('hi')
        foo.connect()
        self.assertEqual(foo._ftp, 'hi')

    def test_disconnect(self):
        # test disconnect where _ftp is None
        foo = FtpFileTransfer(None)
        foo.disconnect()

        # test disconnect where _ftp is set and
        # so is _alt_ftp_con
        foo = FtpFileTransfer(None)
        foo.set_connection('hi')
        foo.connect()
        foo.disconnect()

        # test disconnect where _ftp is set and
        # _alt_ftp_con is None
        foo = FtpFileTransfer(None)
        mockftp = MockFtp()
        mockftp.close = Mock(side_effect=Exception())
        foo.set_connection(mockftp)
        foo.connect()
        foo._alt_ftp_con = None
        foo.disconnect()
        mockftp.close.assert_any_call()

    def test_upload_file_direct_none_passed_as_file(self):
        foo = FtpFileTransfer(None)

        self.assertEqual(foo.upload_file_direct(None, '/remote_dir',
                                                'remote_filename'), False)

        self.assertEqual(foo.get_error_msg(), 'File passed in is None')

    def test_upload_file_direct_file_is_not_a_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileTransfer(None)
            noexist = os.path.join(temp_dir, 'noexist')
            self.assertEqual(foo.upload_file_direct(noexist,
                                                    '/remote_dir',
                                                    'remote_filename'), False)

            self.assertEqual(foo.get_error_msg(), noexist +
                             ' is not a file')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_direct_remote_dir_is_none(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = FtpFileTransfer(None)
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
            foo = FtpFileTransfer(None)
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
            foo = FtpFileTransfer(None)
            foo.set_connection(mockftp)
            foo.connect()
            afile = os.path.join(temp_dir, 'afile')
            open(afile, 'a').close()
            self.assertEqual(foo.upload_file_direct(afile, '/foo', 'name'),
                             False)
            self.assertEqual(foo.get_error_msg(), 'Unable to upload ' + afile +
                             ' to /foo/name : hi')
            foo.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_direct_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockftp = MockFtp()
            mockftp.put = Mock(return_value=3)
            foo = FtpFileTransfer(None)
            foo.set_connection(mockftp)
            foo.connect()
            afile = os.path.join(temp_dir, 'afile')
            open(afile, 'a').close()
            self.assertEqual(foo.upload_file_direct(afile, '/foo', 'name'),
                             True)
            self.assertEqual(foo.get_error_msg(), None)
            self.assertEqual(foo.get_upload_summary(), '1 (3 bytes) files '
                                                       'uploaded in 0 seconds'
                                                       ' to host Unset:')
            foo.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_where_file_is_none(self):

        # test where file is None
        ftpspec = ["put"]
        mockftp = Mock(spec=ftpspec)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        foo._upload_file(None)

    def test_upload_file_where_file_does_not_exist(self):
        temp_dir = tempfile.mkdtemp()
        try:
            ftpspec = ["put"]
            mockftp = Mock(spec=ftpspec)
            foo = FtpFileTransfer(None)
            foo.set_connection(mockftp)
            foo.connect()
            foo._upload_file(os.path.join(temp_dir, 'nonexist'))
            self.assertEqual(mockftp.call_count, 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_upload_file_that_raises_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:

            mockftp = MockFtp()
            mockftp.put = Mock(side_effect=IOError('hi'))
            foo = FtpFileTransfer(None)
            foo.set_connection(mockftp)
            foo.connect()
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
            foo = FtpFileTransfer(None)
            foo.set_connection(mockftp)
            foo.set_remote_dir('/remote')
            foo.connect()
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
        foo = FtpFileTransfer(None)
        self.assertEqual(foo.upload_files(None), True)
        self.assertEqual(foo.get_error_msg(),
                         'List of files passed in was None')
        self.assertEqual(foo.get_upload_summary(),
                         'List of files passed in was None\n'
                         '0 (0 bytes) files uploaded in 0 '
                         'seconds to host Unset:')

    def test_upload_files_file_list_is_empty(self):
        foo = FtpFileTransfer(None)
        self.assertEqual(foo.upload_files([]), True)
        self.assertEqual(foo.get_error_msg(),
                         'No files to upload')
        self.assertEqual(foo.get_upload_summary(),
                         'No files to upload\n'
                         '0 (0 bytes) files uploaded in 0 '
                         'seconds to host Unset:')

    def test_upload_files_one_valid_file(self):
        mockftp = MockFtp()
        mockftp.put = Mock(return_value=3)
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)

        foo.set_remote_dir('/remote')
        foo.set_host('hosty')
        temp_dir = tempfile.mkdtemp()
        try:
            foo.connect()
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
            foo.disconnect()
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
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)

        foo.set_remote_dir('/remote')
        foo.set_host('hosty')
        temp_dir = tempfile.mkdtemp()
        try:
            foo.connect()
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
            foo.disconnect()
            self.assertEqual(mockftp.put.call_count, 2)

            mockftp.close.assert_not_called()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_files_where_file_upload_raises_exception(self):
        mockftp = MockFtp()
        mockftp.put = Mock(side_effect=IOError('hi'))
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.set_remote_dir('/remote')
        foo.set_host('hosty')
        temp_dir = tempfile.mkdtemp()
        try:
            foo.connect()
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
            foo.disconnect()
            mockftp.put.assert_called_with(valid_file,
                                           os.path.normpath('/remote' +
                                                            valid_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_get_upload_summary_remotedir_not_set(self):
        foo = FtpFileTransfer(None)
        sum = foo.get_upload_summary()
        self.assertEqual(sum, '0 (0 bytes) files uploaded in 0 seconds to'
                              ' host Unset:')

    def test_delete_file_none_file(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.delete_file(None))
        self.assertEqual(foo.get_error_msg(), 'remote_file None')

    def test_delete_file_not_connected(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.delete_file('/./somefile'))
        self.assertEqual(foo.get_error_msg(), "Unable to delete /somefile : "
                                              "'NoneType' object has no "
                                              "attribute 'delete'")

    def test_delete_file_success(self):
        mockftp = MockFtp()
        mockftp.delete = Mock(return_value='hello')
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        self.assertTrue(foo.delete_file('/a/b'))
        self.assertEqual(foo.get_error_msg(), None)
        foo.disconnect()
        mockftp.delete.assert_called_with('/a/b')

    def test_delete_file_fail(self):
        mockftp = MockFtp()
        mockftp.delete = Mock(side_effect=IOError('error'))
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        self.assertFalse(foo.delete_file('/a/b'))
        self.assertEqual(foo.get_error_msg(), 'Unable to delete /a/b : error')
        foo.disconnect()
        mockftp.delete.assert_called_with('/a/b')

    def test_download_file_remote_file_params_none(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.download_file(None, None))
        self.assertEqual(foo.get_error_msg(), 'remote_file None')

        self.assertFalse(foo.download_file('/hi', None))
        self.assertEqual(foo.get_error_msg(), 'local_file None')

    def test_download_file_not_connected(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.download_file('/remote/bye', '/local/hi'))
        self.assertEqual(foo.get_error_msg(), "Unable to download /remote/bye "
                                              "to /local/hi : 'NoneType' "
                                              "object has no attribute 'get'")

    def test_download_file_success(self):
        mockftp = MockFtp()
        mockftp.get = Mock()
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        self.assertTrue(foo.download_file('/a/b', '/c/d'))
        self.assertEqual(foo.get_error_msg(), None)
        foo.disconnect()
        mockftp.get.assert_called_with('/a/b', local='/c/d')

    def test_download_file_fail(self):
        mockftp = MockFtp()
        mockftp.get = Mock(side_effect=IOError('error'))
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        self.assertFalse(foo.download_file('/a/b', '/c/d'))
        self.assertEqual(foo.get_error_msg(), 'Unable to download /a/b to '
                                              '/c/d : error')
        foo.disconnect()
        mockftp.get.assert_called_with('/a/b', local='/c/d')

    def test_list_dirs_remote_dir_none(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.list_dirs(None))
        self.assertEqual(foo.get_error_msg(), 'remote_dir None')

    def test_list_dirs_remote_dir_not_connected(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.list_dirs('/foo'))
        self.assertEqual(foo.get_error_msg(), "Unable to get directory list "
                                              "for /foo : 'NoneType' object "
                                              "has no attribute 'list'")

    def test_list_dirs_success(self):
        mockftp = MockFtp()
        mockftp.list = Mock(return_value=[{'directory': 'd', 'name': '.'},
                                          {'directory': 'd', 'name': '..'},
                                          {'directory': 'd', 'name': 'foo'},
                                          {'directory': '-', 'name': 'file'}])
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        dirlist = foo.list_dirs('/foo')
        self.assertTrue(len(dirlist) == 1)
        self.assertEqual(dirlist[0], 'foo')
        foo.disconnect()
        mockftp.list.assert_called_with('/foo', extra=True)

    def test_list_dirs_fail(self):
        mockftp = MockFtp()
        mockftp.list = Mock(side_effect=IOError('error'))
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        dirlist = foo.list_dirs('/foo2')
        self.assertTrue(dirlist is None)
        foo.disconnect()
        self.assertEqual(foo.get_error_msg(), 'Unable to get directory list '
                                              'for /foo2 : error')
        mockftp.list.assert_called_with('/foo2', extra=True)

    def test_list_files_remote_dir_none(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.list_files(None))
        self.assertEqual(foo.get_error_msg(), 'remote_dir None')

    def test_list_files_remote_dir_not_connected(self):
        foo = FtpFileTransfer(None)
        self.assertFalse(foo.list_files('/foo'))
        self.assertEqual(foo.get_error_msg(), "Unable to get file list "
                                              "for /foo : 'NoneType' object "
                                              "has no attribute 'list'")

    def test_list_files_success(self):
        mockftp = MockFtp()
        mockftp.list = Mock(return_value=[{'directory': 'd', 'name': '.'},
                                          {'directory': 'd', 'name': '..'},
                                          {'directory': 'd', 'name': 'foo'},
                                          {'directory': '-', 'name': 'file'}])
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        filelist = foo.list_files('/foo')
        self.assertTrue(len(filelist) == 1)
        self.assertEqual(filelist[0], 'file')
        foo.disconnect()
        mockftp.list.assert_called_with('/foo', extra=True)

    def test_list_files_fail(self):
        mockftp = MockFtp()
        mockftp.list = Mock(side_effect=IOError('error'))
        mockftp.close = Mock(return_value=None)
        foo = FtpFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        filelist = foo.list_files('/foo2')
        self.assertTrue(filelist is None)
        foo.disconnect()
        self.assertEqual(foo.get_error_msg(), 'Unable to get file list '
                                              'for /foo2 : error')
        mockftp.list.assert_called_with('/foo2', extra=True)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
