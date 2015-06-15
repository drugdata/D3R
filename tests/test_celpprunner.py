#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_celpprunner
----------------------------------

Tests for `celpprunner` module.
"""

import unittest
import tempfile
import argparse
import logging
import os.path

from argparse import ArgumentError

from d3r import celpprunner
from d3r.task import D3RParameters

class TestCelppRunner(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_lock(self):
        tempDir = tempfile.mkdtemp()
        theargs = D3RParameters()
        theargs.celppdir=tempDir
        theargs.stage='blast'

        # get the lock file which should work
        lock = celpprunner._get_lock(theargs)
        expectedLockFile = os.path.join(tempDir,'celpprunner.blast.lockpid')
        self.assertTrue(os.path.isfile(expectedLockFile))

        # try getting lock again which should also work
        lock = celpprunner._get_lock(theargs)

        lock.release()
        self.assertFalse(os.path.isfile(expectedLockFile))
        os.rmdir(tempDir)       
         
    def test_setup_logging(self):
        logger = logging.getLogger('funlogger')
        theargs = D3RParameters()
        theargs.logLevel = 'DEBUG'
        celpprunner._setup_logging(theargs)
        self.assertEqual(logger.getEffectiveLevel(),30)
       

    def test_parse_arguments(self):
        theargs = ['--stage','blast','foo']
        
        result = celpprunner._parse_arguments('hi',theargs)
        self.assertEqual(result.stage,'blast')
        self.assertEqual(result.celppdir,'foo')
        

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
