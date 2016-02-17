#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import os.path

"""
test_uploader
----------------------------------

Tests for `uploader` module.
"""


class TestFtpFileUploader(unittest.TestCase):

    def setUp(self):
        pass


    def test_constructor(self):
        self.assertEqual(1, 2)

        # test passing None

        # test passing empty config file

        # test passing invalid config file

        # test passing partial config file

        # test passing valid config file

    # test getters and setters

    # test connect

    # test disconnect

    # test _upload_file

    # test upload_files

    # test get_upload_summary

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
