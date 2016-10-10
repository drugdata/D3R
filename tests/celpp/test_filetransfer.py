#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import shutil
import os.path

from d3r.celpp.filetransfer import FileTransfer
from d3r.celpp.filetransfer import InvalidFtpConfigException
"""
test_uploader
----------------------------------

Tests for `uploader` module.
"""


class MockFtp(object):
    """Dummy ftp connection
    """


class TestFileTransfer(unittest.TestCase):

    def setUp(self):
        pass

    def test_filetransfer_parse_config(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test _parse_config with None
            foo = FileTransfer(None)
            try:
                foo._parse_config(None)
                self.fail('Expected InvalidFtpConfigException')
            except InvalidFtpConfigException:
                pass

            # test passing invalid config file
            f = open(os.path.join(temp_dir, 'invalid'), 'a')
            f.write('hello\nhow\nare')
            f.flush()
            f.close()
            foo = FileTransfer(os.path.join(temp_dir, 'invalid'))
            self.assertEqual(foo.get_host(), None)
            self.assertEqual(foo.get_password(), None)
            self.assertEqual(foo.get_user(), None)
            self.assertEqual(foo.get_remote_dir(), '')
            self.assertEqual(foo.get_remote_challenge_dir(), None)
            self.assertEqual(foo.get_remote_submission_dir(), None)

            # test passing partial config file
            f = open(os.path.join(temp_dir, 'partial'), 'a')
            f.write('#hello\nuser bob@bob.com\nare')
            f.flush()
            f.close()
            foo = FileTransfer(os.path.join(temp_dir, 'partial'))
            self.assertEqual(foo.get_host(), None)
            self.assertEqual(foo.get_password(), None)
            self.assertEqual(foo.get_user(), 'bob@bob.com')
            self.assertEqual(foo.get_remote_dir(), '')
            self.assertEqual(foo.get_remote_challenge_dir(), None)
            self.assertEqual(foo.get_remote_submission_dir(), None)

            # test passing valid config file
            f = open(os.path.join(temp_dir, 'valid'), 'a')
            f.write('host ftp.box.com\nuser bob@bob.com\npass 222\npath /foo')
            f.flush()
            f.close()
            foo = FileTransfer(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_host(), 'ftp.box.com')
            self.assertEqual(foo.get_password(), '222')
            self.assertEqual(foo.get_user(), 'bob@bob.com')
            self.assertEqual(foo.get_remote_dir(), '/foo')
            self.assertEqual(foo.get_remote_challenge_dir(), None)
            self.assertEqual(foo.get_remote_submission_dir(), None)

            # test passing valid config file with challengepath
            f = open(os.path.join(temp_dir, 'valid'), 'a')
            f.write('host ftp.box.com\nuser bob@bob.com\npass 222\npath /foo\n'
                    'challengepath /chall\n')
            f.flush()
            f.close()
            foo = FileTransfer(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_host(), 'ftp.box.com')
            self.assertEqual(foo.get_password(), '222')
            self.assertEqual(foo.get_user(), 'bob@bob.com')
            self.assertEqual(foo.get_remote_dir(), '/foo')
            self.assertEqual(foo.get_remote_challenge_dir(), '/chall')
            self.assertEqual(foo.get_remote_submission_dir(), None)

            # test passing valid config file with challengepath and
            # submissionpath
            f = open(os.path.join(temp_dir, 'valid'), 'a')
            f.write('host ftp.box.com\nuser bob@bob.com\npass 222\npath /foo\n'
                    'challengepath /chall\nsubmissionpath /submit\n')
            f.flush()
            f.close()
            foo = FileTransfer(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_host(), 'ftp.box.com')
            self.assertEqual(foo.get_password(), '222')
            self.assertEqual(foo.get_user(), 'bob@bob.com')
            self.assertEqual(foo.get_remote_dir(), '/foo')
            self.assertEqual(foo.get_remote_challenge_dir(), '/chall')
            self.assertEqual(foo.get_remote_submission_dir(), '/submit')

        finally:
            shutil.rmtree(temp_dir)

    def test_constructor(self):
        # test passing None which means we set everything manually
        foo = FileTransfer(None)

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
                foo = FileTransfer(os.path.join(temp_dir, 'doesnotexist'))
                self.fail('Expected IOError')
            except IOError:
                pass

            # test passing empty config file
            open(os.path.join(temp_dir, 'empty'), 'a').close()
            foo = FileTransfer(os.path.join(temp_dir, 'empty'))
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
        foo = FileTransfer(None)
        self.assertTrue(foo.connect())

    def test_disconnect(self):
        foo = FileTransfer(None)
        foo.disconnect()

    def test_ftp_upload_file_direct(self):
        foo = FileTransfer(None)

        self.assertEqual(foo.upload_file_direct(None, '/remote_dir',
                                                'remote_filename'), False)

        self.assertEqual(foo.get_error_msg(),
                         'upload_file_direct not implemented')

    def test_upload_files(self):
        foo = FileTransfer(None)
        self.assertEqual(foo.upload_files(None), False)
        self.assertEqual(foo.get_error_msg(),
                         'upload_files not implemented')
        self.assertEqual(foo.get_upload_summary(),
                         'upload_files not implemented\n0 (0 bytes) files '
                         'uploaded in 0 seconds to host Unset:')

    def test_get_upload_summary(self):
        foo = FileTransfer(None)
        sum = foo.get_upload_summary()
        self.assertEqual(sum, '0 (0 bytes) files uploaded in 0 seconds to'
                              ' host Unset:')

    def test_delete_file_none_file(self):
        foo = FileTransfer(None)
        self.assertFalse(foo.delete_file(None))
        self.assertEqual(foo.get_error_msg(), 'delete_file not implemented')

    def test_download_file(self):
        foo = FileTransfer(None)
        self.assertFalse(foo.download_file(None, None))
        self.assertEqual(foo.get_error_msg(), 'download_file not implemented')

    def test_list_dirs(self):
        foo = FileTransfer(None)
        self.assertFalse(foo.list_dirs(None))
        self.assertEqual(foo.get_error_msg(), 'list_dirs not implemented')

    def test_list_files(self):
        foo = FileTransfer(None)
        self.assertFalse(foo.list_files(None))
        self.assertEqual(foo.get_error_msg(), 'list_files not implemented')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
