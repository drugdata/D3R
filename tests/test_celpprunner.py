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
from datetime import date

from d3r import celpprunner
from d3r.celpp.task import D3RParameters
from d3r.celpp import util
from d3r.celpp.task import D3RTask


class DummyTask(D3RTask):
    """Dummy Task used for tests below
    """
    def __init__(self, path, args, error_message, can_run_return,
                 can_run_exception, run_exception):
        """Sets up a dummy task
        """
        super(DummyTask, self).__init__(path, args)
        self.set_error(error_message)
        self._can_run_return = can_run_return
        self._can_run_exception = can_run_exception
        self._run_exception = run_exception
        self._run_count = 0

    def can_run(self):
        if self._can_run_exception is not None:
            raise self._can_run_exception

        return self._can_run_return

    def run(self):
        self._run_count += 1
        if self._run_exception is not None:
            raise self._run_exception


class TestCelppRunner(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_lock(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir

            # get the lock file which should work
            lock = celpprunner._get_lock(theargs, 'blast')
            expectedLockFile = os.path.join(temp_dir,
                                            'celpprunner.blast.lockpid')
            self.assertTrue(os.path.isfile(expectedLockFile))

            # try getting lock again which should also work
            lock = celpprunner._get_lock(theargs, 'blast')

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
        logger.debug('test')

        theargs.loglevel = 'DEBUG'
        celpprunner._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.task').getEffectiveLevel(),
                         logging.DEBUG)
        self.assertEqual(theargs.numericloglevel, logging.DEBUG)

        theargs.loglevel = 'WARNING'
        celpprunner._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.task').getEffectiveLevel(),
                         logging.WARNING)
        self.assertEqual(theargs.numericloglevel, logging.WARNING)

        theargs.loglevel = 'ERROR'
        celpprunner._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.task').getEffectiveLevel(),
                         logging.ERROR)
        self.assertEqual(theargs.numericloglevel, logging.ERROR)

        theargs.loglevel = 'CRITICAL'
        celpprunner._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.task').getEffectiveLevel(),
                         logging.CRITICAL)
        self.assertEqual(theargs.numericloglevel, logging.CRITICAL)

    def test_parse_arguments(self):
        theargs = ['--stage', 'blast', 'foo']
        result = celpprunner._parse_arguments('hi', theargs)
        self.assertEqual(result.stage, 'blast')
        self.assertEqual(result.celppdir, 'foo')
        self.assertEqual(result.blastdir, None)
        self.assertEqual(result.email, None)
        self.assertEqual(result.loglevel, 'WARNING')
        self.assertEqual(result.blastnfilter, 'blastnfilter.py')
        self.assertEqual(result.pdbprep, 'pdbprep.py')

        theargs = ['foo', '--stage', 'dock', '--email', 'b@b.com,h@h',
                   '--blastdir', 'b', '--log', 'ERROR',
                   '--blastnfilter', '/bin/blastnfilter.py',
                   '--pdbprep', '/bin/pdbprep.py',
                   '--postanalysis', '/bin/postanalysis.py']
        result = celpprunner._parse_arguments('hi', theargs)
        self.assertEqual(result.stage, 'dock')
        self.assertEqual(result.celppdir, 'foo')
        self.assertEqual(result.blastdir, 'b')
        self.assertEqual(result.email, 'b@b.com,h@h')
        self.assertEqual(result.loglevel, 'ERROR')
        self.assertEqual(result.blastnfilter, '/bin/blastnfilter.py')
        self.assertEqual(result.pdbprep, '/bin/pdbprep.py')
        self.assertEquals(result.postanalysis, '/bin/postanalysis.py')

    def test_run_tasks_passing_none_and_empty_list(self):
        self.assertEquals(celpprunner.run_tasks(None), 3)
        task_list = []
        self.assertEquals(celpprunner.run_tasks(task_list), 2)

    def test_run_one_successful_task(self):
        success_task = DummyTask(D3RParameters(), 'foo', None, True, None,
                                 None)
        success_task.set_name('dummy')
        task_list = []
        task_list.append(success_task)
        self.assertEquals(celpprunner.run_tasks(task_list), 0)

    def test_run_one_fail_task_with_error_message(self):
        task = DummyTask(D3RParameters(), 'foo', 'someerror', True, None, None)
        task.set_name('dummy')
        task_list = []
        task_list.append(task)
        self.assertEquals(celpprunner.run_tasks(task_list), 1)
        self.assertEquals(task.get_error(), 'someerror')

    def test_run_one_fail_task_with_exception_and_no_message(self):
        task = DummyTask(D3RParameters(), 'foo', None, True,
                         None, Exception('hi'))
        task.set_name('dummy')
        task_list = []
        task_list.append(task)
        self.assertEquals(celpprunner.run_tasks(task_list), 1)
        self.assertEquals(task.get_error(),
                          'Caught Exception running task: hi')

    def test_run_two_tasks_success(self):
        task_list = []
        task = DummyTask(D3RParameters(), 'foo', None, True, None, None)
        task.set_name('dummy')
        task_list.append(task)
        task_list.append(task)

        self.assertEquals(celpprunner.run_tasks(task_list), 0)
        self.assertEquals(task._run_count, 2)

    def test_run_two_tasks_second_task_has_error(self):
        task_list = []
        task = DummyTask(D3RParameters(), 'foo', None, True, None, None)
        task.set_name('dummy')
        task_list.append(task)

        task_two = DummyTask(D3RParameters(), 'foo', None, True,
                             None, Exception('hi'))
        task_two.set_name('dummy')
        task_list.append(task_two)

        self.assertEquals(celpprunner.run_tasks(task_list), 1)
        self.assertEquals(task._run_count, 1)
        self.assertEquals(task_two._run_count, 1)
        self.assertEquals(task_two.get_error(),
                          'Caught Exception running task: hi')

    def test_run_two_tasks_first_task_has_error(self):
        task_list = []
        task = DummyTask(D3RParameters(), 'foo', None, True, None,
                         Exception('hi'))
        task.set_name('dummy')
        task_list.append(task)

        task_two = DummyTask(D3RParameters(), 'foo', None, True, None,
                             None)
        task_two.set_name('dummy')
        task_list.append(task_two)

        self.assertEquals(celpprunner.run_tasks(task_list), 1)
        self.assertEquals(task.get_error(),
                          'Caught Exception running task: hi')

        self.assertEquals(task._run_count, 1)
        self.assertEquals(task_two._run_count, 0)

    def test_get_task_list_for_stage_with_invalid_stage_name(self):

        try:
            celpprunner.get_task_list_for_stage(D3RParameters(), None)
            self.fail('Expected exception')
        except NotImplementedError as e:
            self.assertEquals(e.message, 'stage_name is None')

        try:
            celpprunner.get_task_list_for_stage(D3RParameters(), '')
            self.fail('Expected exception')
        except NotImplementedError as e:
            self.assertEquals(e.message, 'uh oh no tasks for  stage')

        try:
            celpprunner.get_task_list_for_stage(D3RParameters(), 'foo')
            self.fail('Expected exception')
        except NotImplementedError as e:
            self.assertEquals(e.message, 'uh oh no tasks for foo stage')

    def test_get_task_list_for_stage_with_valid_stages(self):
        params = D3RParameters()
        params.latest_weekly = 'foo'
        task_list = celpprunner.get_task_list_for_stage(params, 'blast')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', 'stage.2.blastnfilter'))

        task_list = celpprunner.get_task_list_for_stage(params, 'pdbprep')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', 'stage.3.pdbprep'))

        task_list = celpprunner.get_task_list_for_stage(params, 'import')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', 'stage.1.dataimport'))

    def test_run_stages_no_weekly_datasetfound(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir
            self.assertEqual(celpprunner.run_stages(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_invalid_stage(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir
            theargs.stage = 'foo'
            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))
            try:
                celpprunner.run_stages(theargs)
            except NotImplementedError as e:
                self.assertEquals(e.message, 'uh oh no tasks for foo stage')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_blast_stage_data_import_missing(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir
            theargs.stage = 'blast'
            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'complete'), 'a').close()
            theargs.blastdir = temp_dir

            self.assertEquals(celpprunner.run_stages(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_blast(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            theargs.stage = 'blast'
            theargs.pdbdb = '/pdbdb'
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'complete'), 'a').close()
            theargs.blastdir = temp_dir
            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        'stage.1.dataimport')
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir, 'complete'), 'a').close()

            theargs.blastnfilter = 'echo'
            theargs.postanalysis = 'true'
            self.assertEqual(celpprunner.run_stages(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_blast_has_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            theargs.stage = 'blast'
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'error'), 'a').close()
            theargs.blastdir = temp_dir
            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))
            self.assertEqual(celpprunner.run_stages(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_blast_and_pdbprep_no_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.pdbdb = '/pdbdb'
            theargs.celppdir = os.path.join(temp_dir)
            theargs.stage = 'blast,pdbprep'
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'complete'), 'a').close()
            theargs.blastdir = temp_dir
            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        'stage.1.dataimport')
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir, 'complete'), 'a').close()

            theargs.blastnfilter = 'echo'
            theargs.postanalysis = 'true'
            theargs.pdbprep = 'echo'
            self.assertEqual(celpprunner.run_stages(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_blast_and_pdbprep_blast_has_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            theargs.stage = 'blast,pdbprep'
            os.mkdir(os.path.join(temp_dir, 'current'))
            open(os.path.join(temp_dir, 'current', 'complete'), 'a').close()
            theargs.blastdir = temp_dir
            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        'stage.1.dataimport')
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir, 'error'), 'a').close()
            theargs.blastnfilter = 'echo'
            theargs.pdbprep = 'echo'
            self.assertEqual(celpprunner.run_stages(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_createweekdir_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            theargs.createweekdir = True
            theargs.stage = ''
            d = date.today()
            celp_week = util.get_celpp_week_of_year_from_date(d)
            try:
                self.assertEquals(celpprunner.run_stages(theargs), 0)
                self.fail('Expected NotImplementedError')
            except NotImplementedError:
                pass

            expected_dir = os.path.join(temp_dir, str(celp_week[1]),
                                        'dataset.week.' +
                                        str(celp_week[0]))
            self.assertEquals(os.path.isdir(expected_dir), True)

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
