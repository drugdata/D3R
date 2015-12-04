__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_proteinligprep
--------------------------------

Tests for `proteinligprep` module.
"""

import shutil
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.scoring import ScoringTaskFactory
from d3r.celpp.scoring import ScoringTask
from d3r.celpp.scoring import PathNotDirectoryError


class TestScoring(unittest.TestCase):
    def setUp(self):
        pass

    def test_scoringtaskfactory_constructor(self):
        params = D3RParameters()
        params.hi = True
        stf = ScoringTaskFactory('/foo', params)
        self.assertEquals(stf.get_args().hi, True)
        self.assertEquals(stf.get_path(), '/foo')

    def test_get_scoring_tasks_invalid_path(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            path = os.path.join(temp_dir, 'doesnotexist')
            stf = ScoringTaskFactory(path, params)
            try:
                task_list = stf.get_scoring_tasks()
                self.fail('Expected PathNotDirectoryError')
            except PathNotDirectoryError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_get_scoring_tasks_on_empty_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            stf = ScoringTaskFactory(temp_dir, params)
            task_list = stf.get_scoring_tasks()
            self.assertEquals(len(task_list), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_scoringtask_constructor(self):
        params = D3RParameters()
        # no dock task found so it cannot run
        docktask = D3RTask('/blah', params)
        docktask.set_name('foo')
        docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)

        scoring = ScoringTask('/blah', 'foo.scoring',
                              docktask, params)
        self.assertEquals(scoring.get_name(),'foo.scoring')
        self.assertEquals(scoring.get_stage(), 5)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # no dock task found so it cannot run
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)

            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            self.assertEqual(scoring.can_run(), False)
            self.assertEqual(scoring.get_error(),
                             'foo task has notfound status')

            # docktask  running
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.START_FILE),
                 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            self.assertEqual(scoring.can_run(), False)
            self.assertEqual(scoring.get_error(),
                             'foo task has start status')

            # docktask failed
            error_file = os.path.join(docktask.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            self.assertEqual(scoring.can_run(), False)
            self.assertEqual(scoring.get_error(),
                             'foo task has error status')

            # docktask success
            os.remove(error_file)
            open(os.path.join(docktask.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            self.assertEqual(scoring.can_run(), True)
            self.assertEqual(scoring.get_error(), None)

            # scoring task exists already
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.create_dir()
            self.assertEqual(scoring.can_run(), False)
            self.assertEqual(scoring.get_error(),
                             scoring.get_dir_name() +
                             ' already exists and status is unknown')

            # scoring task already complete
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            open(os.path.join(scoring.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(scoring.can_run(), False)
            self.assertEqual(scoring.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # return immediately cause can_run is false
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.run()
            self.assertEqual(scoring.get_error(),
                             'foo task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_scoring_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.run()
            self.assertEqual(scoring.get_error(),
                             'scoring not set')
            # test files get created
            self.assertEqual(os.path.isdir(scoring.get_dir()),
                             True)
            errfile = os.path.join(scoring.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_pdbdb_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.scoring = 'true'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.run()
            self.assertEqual(scoring.get_error(),
                             'pdbdb not set')
            # test files get created
            self.assertEqual(os.path.isdir(scoring.get_dir()),
                             True)
            errfile = os.path.join(scoring.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)
    def test_run_fails_cause_scoring_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.scoring = 'false'
            params.pdbdb = '/data/pdb'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.run()
            self.assertEqual(scoring.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(scoring.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(scoring.get_dir(),
                                   'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(scoring.get_dir(),
                                   'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_scoring_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.scoring = '/bin/doesnotexist'
            params.pdbdb = '/data/pdb'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.run()
            self.assertEqual(scoring.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --pdbdb /data/pdb ' +
                             '--dockdir ' +
                             docktask.get_dir() + ' --outdir ' +
                             scoring.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(scoring.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.scoring = 'true'
            params.pdbdb = '/data/pdb'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            scoring = ScoringTask(temp_dir, 'foo.scoring',
                                  docktask, params)
            scoring.run()
            self.assertEqual(scoring.get_error(), None)
            # test files get created
            errfile = os.path.join(scoring.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(scoring.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(scoring.get_dir(),
                                   'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(scoring.get_dir(),
                                   'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
