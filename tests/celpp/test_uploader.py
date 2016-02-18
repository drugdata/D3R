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

            # test getters and setters
            foo.set_connect_timeout(10)
            foo.set_ftp_connection('hi')
            foo.set_ftp_host('host')
            foo.set_ftp_password('pass')
            foo.set_ftp_remote_dir('/remote')
            foo.set_ftp_user('user')

            self.assertEqual(foo.get_connect_timeout(), 10)
            self.assertEqual(foo.get_ftp_host(), 'host')
            self.assertEqual(foo.get_ftp_password(), 'pass')
            self.assertEqual(foo.get_ftp_remote_dir(), '/remote')
            self.assertEqual(foo.get_ftp_user(), 'user')
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
            f.write('host ftp.box.com\nuser bob@bob.com\npass 222\npath /foo')
            f.flush()
            f.close()
            foo = FtpFileUploader(os.path.join(temp_dir, 'valid'))
            self.assertEqual(foo.get_ftp_host(), 'ftp.box.com')
            self.assertEqual(foo.get_ftp_password(), '222')
            self.assertEqual(foo.get_ftp_user(), 'bob@bob.com')
            self.assertEqual(foo.get_ftp_remote_dir(), '/foo')

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
        foo.set_ftp_connection('hi')
        foo._connect()
        foo._alt_ftp_con = None
        foo._disconnect()



    def test_upload_file(self):
        self.assertEqual(1, 2)



            # test upload_files

            # test get_upload_summary

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
