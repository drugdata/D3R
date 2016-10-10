#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import shutil
import os.path
from mock import Mock

from d3r.celpp.filetransfer import WebDavFileTransfer

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
        foo = WebDavFileTransfer(None)

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
                foo = WebDavFileTransfer(os.path.join(temp_dir,
                                                      'doesnotexist'))
                self.fail('Expected IOError')
            except IOError:
                pass

            # test passing empty config file
            open(os.path.join(temp_dir, 'empty'), 'a').close()
            foo = WebDavFileTransfer(os.path.join(temp_dir, 'empty'))
            self.assertEqual(foo.get_host(), None)
            self.assertEqual(foo.get_password(), None)
            self.assertEqual(foo.get_user(), None)
            self.assertEqual(foo.get_remote_dir(), '')
            self.assertEqual(foo.get_remote_challenge_dir(), None)

            # test getters and setters
            foo.set_connect_timeout(10)
            foo.set_connection('hi')
            foo.set_host('host')
            foo.set_password('pass')
            foo.set_remote_dir('/remote')
            foo.set_user('user')
            foo.set_remote_challenge_dir('/chall')

            self.assertEqual(foo.get_connect_timeout(), 10)
            self.assertEqual(foo.get_host(), 'host')
            self.assertEqual(foo.get_password(), 'pass')
            self.assertEqual(foo.get_remote_dir(), '/remote')
            self.assertEqual(foo.get_user(), 'user')
            self.assertEqual(foo.get_remote_challenge_dir(), '/chall')
        finally:
            shutil.rmtree(temp_dir)

    def test_connect(self):
        # test where alt_ftp_con is None and we use easywebdav
        foo = WebDavFileTransfer(None)
        foo.set_host('doesnotexist')
        foo.set_user('user')
        foo.set_password('')
        foo.set_remote_dir('/remote')
        foo.set_connect_timeout(0)
        self.assertTrue(foo.connect())

        # test where alt_ftp_con is set
        foo = WebDavFileTransfer(None)
        foo.set_connection('hi')
        foo.connect()
        self.assertEqual(foo._ftp, 'hi')

    def test_disconnect(self):
        # test disconnect where _ftp is None
        foo = WebDavFileTransfer(None)
        foo.disconnect()

        # test disconnect where _ftp is set and
        # so is _alt_ftp_con
        foo = WebDavFileTransfer(None)
        foo.set_connection('hi')
        foo.connect()
        foo.disconnect()

        # test disconnect where _ftp is set and
        # _alt_ftp_con is None
        foo = WebDavFileTransfer(None)
        foo.connect()
        foo._alt_ftp_con = None
        foo.disconnect()

    def test_upload_file_direct_none_passed_as_file(self):
        foo = WebDavFileTransfer(None)

        self.assertEqual(foo.upload_file_direct(None, '/remote_dir',
                                                'remote_filename'), False)

        self.assertEqual(foo.get_error_msg(), 'File passed in is None')

    def test_upload_file_direct_file_is_not_a_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            foo = WebDavFileTransfer(None)
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
            foo = WebDavFileTransfer(None)
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
            foo = WebDavFileTransfer(None)
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
            mockftp.upload = Mock(side_effect=IOError('hi'))
            foo = WebDavFileTransfer(None)
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
            mockftp.upload = Mock()
            foo = WebDavFileTransfer(None)
            foo.set_connection(mockftp)
            foo.connect()
            fakefile = os.path.join(temp_dir, 'foo')
            f = open(fakefile, 'w')
            f.write('123')
            f.flush()
            f.close()
            self.assertEqual(foo.upload_file_direct(fakefile, '/foo', 'name'),
                             True)
            self.assertEqual(foo.get_error_msg(), None)
            self.assertEqual(foo.get_upload_summary(), '1 (3 bytes) files '
                                                       'uploaded in 0 seconds'
                                                       ' to host Unset:')
            foo.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def test_get_upload_summary_remotedir_not_set(self):
        foo = WebDavFileTransfer(None)
        sum = foo.get_upload_summary()
        self.assertEqual(sum, '0 (0 bytes) files uploaded in 0 seconds to'
                              ' host Unset:')

    def test_download_file_remote_file_params_none(self):
        foo = WebDavFileTransfer(None)
        self.assertFalse(foo.download_file(None, None))
        self.assertEqual(foo.get_error_msg(), 'remote_file None')

        self.assertFalse(foo.download_file('/hi', None))
        self.assertEqual(foo.get_error_msg(), 'local_file None')

    def test_download_file_not_connected(self):
        foo = WebDavFileTransfer(None)
        self.assertFalse(foo.download_file('/remote/bye', '/local/hi'))
        self.assertEqual(foo.get_error_msg(), "Unable to download /remote/bye "
                                              "to /local/hi : 'NoneType' "
                                              "object has no attribute "
                                              "'download'")

    def test_download_file_success(self):
        mockftp = MockFtp()
        mockftp.download = Mock()
        foo = WebDavFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        self.assertTrue(foo.download_file('/a/b', '/c/d'))
        self.assertEqual(foo.get_error_msg(), None)
        foo.disconnect()
        mockftp.download.assert_called_with('/a/b', '/c/d')

    def test_download_file_fail(self):
        mockftp = MockFtp()
        mockftp.download = Mock(side_effect=IOError('error'))
        foo = WebDavFileTransfer(None)
        foo.set_connection(mockftp)
        foo.connect()
        self.assertFalse(foo.download_file('/a/b', '/c/d'))
        self.assertEqual(foo.get_error_msg(), 'Unable to download /a/b to '
                                              '/c/d : error')
        foo.disconnect()
        mockftp.download.assert_called_with('/a/b', '/c/d')

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
