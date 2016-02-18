#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import shutil
import os.path

from d3r.celpp.uploader import FtpFileUploader
from d3r.celpp.uploader import InvalidFtpConfigException
"""
test_uploader
----------------------------------

Tests for `uploader` module.
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
        self.assertEqual(foo.get_ftp_remote_dir(), None)

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
            self.assertEqual(foo.get_ftp_remote_dir(), None)

            # test passing invalid config file
            f = open(os.path.join(temp_dir, 'invalid'), 'a')
            f.write('hello\nhow\nare')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'invalid'))
            self.assertEqual(foo.get_ftp_host(), None)
            self.assertEqual(foo.get_ftp_password(), None)
            self.assertEqual(foo.get_ftp_user(), None)
            self.assertEqual(foo.get_ftp_remote_dir(), None)

            # test passing partial config file
            f = open(os.path.join(temp_dir, 'partial'), 'a')
            f.write('#hello\nuser bob@bob.com\nare')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'partial'))
            self.assertEqual(foo.get_ftp_host(), None)
            self.assertEqual(foo.get_ftp_password(), None)
            self.assertEqual(foo.get_ftp_user(), 'bob@bob.com')
            self.assertEqual(foo.get_ftp_remote_dir(), None)

            # test passing valid config file
            f = open(os.path.join(temp_dir, 'valid'), 'a')
            f.write('host  ftp.box.com\nuser bob@bob.com\npass 222\npath /foo')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_ftp_host(), 'ftp.box.com')
            self.assertEqual(foo.get_ftp_password(), '222')
            self.assertEqual(foo.get_ftp_user(), 'bob@bob.com')
            self.assertEqual(foo.get_ftp_remote_dir(), '/foo')

            # test getters and setters
            self.assertEqual(1, 2)
            # test connect

            # test disconnect

            # test _upload_file

            # test upload_files

            # test get_upload_summary
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
