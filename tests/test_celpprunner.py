#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_celpprunner
----------------------------------

Tests for `celpprunner` module.
"""

import unittest
import tempfile
import logging
import os
import os.path
import shutil

from d3r import celpprunner
from d3r.task import D3RParameters


class TestCelppRunner(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_lock(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir
            theargs.stage = 'blast'

            # get the lock file which should work
            lock = celpprunner._get_lock(theargs)
            expectedLockFile = os.path.join(temp_dir,
                                            'celpprunner.blast.lockpid')
            self.assertTrue(os.path.isfile(expectedLockFile))

            # try getting lock again which should also work
            lock = celpprunner._get_lock(theargs)

            lock.release()
            self.assertFalse(os.path.isfile(expectedLockFile))
        finally:
            shutil.rmtree(temp_dir)

    def test_setup_logging(self):
        logger = logging.getLogger('funlogger')
        theargs = D3RParameters()
        theargs.loglevel = 'INFO'
        celpprunner._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.task').getEffectiveLevel(),
                         logging.INFO)
        self.assertEqual(theargs.numericloglevel, logging.INFO)

    def test_parse_arguments(self):
        theargs = ['--stage', 'blast', 'foo', '--blastnfilter', 'true']
        result = celpprunner._parse_arguments('hi', theargs)
        self.assertEqual(result.stage, 'blast')
        self.assertEqual(result.celppdir, 'foo')
        self.assertEqual(result.blastdir, None)
        self.assertEqual(result.email, None)
        self.assertEqual(result.loglevel, 'WARNING')
        self.assertEqual(result.blastnfilter, 'true')

        theargs = ['foo', '--stage', 'dock', '--email', 'b@b.com,h@h',
                   '--blastdir', 'b', '--log', 'ERROR',
                   '--blastnfilter', '/bin/blastnfilter.py']
        result = celpprunner._parse_arguments('hi', theargs)
        self.assertEqual(result.stage, 'dock')
        self.assertEqual(result.celppdir, 'foo')
        self.assertEqual(result.blastdir, 'b')
        self.assertEqual(result.email, 'b@b.com,h@h')
        self.assertEqual(result.loglevel, 'ERROR')
        self.assertEqual(result.blastnfilter, '/bin/blastnfilter.py')

    def test_run_stage_no_weekly_datasetfound(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir
            self.assertEqual(celpprunner.run_stage(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stage_dock(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)

            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))

            theargs.stage = 'dock'
            try:
                celpprunner.run_stage(theargs)
                self.fail('Expected NotImplementedError')
            except NotImplementedError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stage_score(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)

            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))

            theargs.stage = 'score'
            try:
                celpprunner.run_stage(theargs)
                self.fail('Expected NotImplementedError')
            except NotImplementedError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stage_blast_data_import_missing(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'complete'), 'a').close()
            theargs.blastdir = temp_dir
            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))

            theargs.stage = 'blast'
            self.assertEqual(celpprunner.run_stage(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stage_blast(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'complete'), 'a').close()
            theargs.blastdir = temp_dir
            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        'stage.1.dataimport')
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir, 'complete'), 'a').close()

            theargs.stage = 'blast'
            theargs.blastnfilter= 'echo'
            self.assertEqual(celpprunner.run_stage(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)


    def test_run_stage_blast_has_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'error'), 'a').close()
            theargs.blastdir = temp_dir
            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))
            theargs.stage = 'blast'
            self.assertEqual(celpprunner.run_stage(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
