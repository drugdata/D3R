__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_evaluation
--------------------------------

Tests for `evaluation` module.
"""
import shutil
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.evaluation import PathNotDirectoryError


class TestEvaluation(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        self.assertEqual('need to ', 'test this')

    def test_evaluationtaskfactory_constructor(self):
        params = D3RParameters()
        params.hi = True
        stf = EvaluationTaskFactory('/foo', params)
        self.assertEquals(stf.get_args().hi, True)
        self.assertEquals(stf.get_path(), '/foo')

    def test_get_evaluation_tasks_invalid_path(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            path = os.path.join(temp_dir, 'doesnotexist')
            stf = EvaluationTaskFactory(path, params)
            try:
                stf.get_evaluation_tasks()
                self.fail('Expected PathNotDirectoryError')
            except PathNotDirectoryError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_empty_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_dir_with_lower_stages_dirs(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            os.mkdir(os.path.join(temp_dir, 'stage.1.dataimport'))
            os.mkdir(os.path.join(temp_dir, 'stage.2.blastnfilter'))

            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_dir_with_webdata_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            os.mkdir(os.path.join(temp_dir,
                                  EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                  EvaluationTaskFactory.WEB_DATA_SUFFIX))

            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_with_valid_algo_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            glidedir = os.path.join(temp_dir,
                                    EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.COMPLETE_FILE), 'a').close()
            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 1)
            self.assertEquals(task_list[0].get_name(), 'glide.evaluation')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_with_two_valid_algo_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            glidedir = os.path.join(temp_dir,
                                    EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.COMPLETE_FILE), 'a').close()

            freddir = os.path.join(temp_dir,
                                   EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                   'fred')
            os.mkdir(freddir)
            open(os.path.join(freddir, D3RTask.COMPLETE_FILE), 'a').close()

            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 2)
            self.assertNotEquals(task_list[0].get_name(),
                                 task_list[1].get_name())

        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_with_one_invalid_algo(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            glidedir = os.path.join(temp_dir,
                                    EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.ERROR_FILE), 'a').close()

            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluation_tasks_on_with_one_valid_and_one_invalid_algo(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            glidedir = os.path.join(temp_dir,
                                    EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.ERROR_FILE), 'a').close()
            freddir = os.path.join(temp_dir,
                                   EvaluationTaskFactory.STAGE_FOUR_PREFIX +
                                   'fred')
            os.mkdir(freddir)
            open(os.path.join(freddir, D3RTask.COMPLETE_FILE), 'a').close()

            stf = EvaluationTaskFactory(temp_dir, params)
            task_list = stf.get_evaluation_tasks()
            self.assertEquals(len(task_list), 1)
            self.assertEquals(task_list[0].get_name(), 'fred.evaluation')
        finally:
            shutil.rmtree(temp_dir)

    def test_evaluationtask_constructor(self):
        params = D3RParameters()
        # no dock task found so it cannot run
        docktask = D3RTask('/blah', params)
        docktask.set_name('foo')
        docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)

        evaluation = EvaluationTask('/blah', 'foo.evaluation',
                                    docktask, params)
        self.assertEquals(evaluation.get_name(), 'foo.evaluation')
        self.assertEquals(evaluation.get_stage(), 5)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # no dock task found so it cannot run
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)

            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            self.assertEqual(evaluation.can_run(), False)
            self.assertEqual(evaluation.get_error(),
                             'foo task has notfound status')

            # docktask  running
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.START_FILE),
                 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            self.assertEqual(evaluation.can_run(), False)
            self.assertEqual(evaluation.get_error(),
                             'foo task has start status')

            # docktask failed
            error_file = os.path.join(docktask.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            self.assertEqual(evaluation.can_run(), False)
            self.assertEqual(evaluation.get_error(),
                             'foo task has error status')

            # docktask success
            os.remove(error_file)
            open(os.path.join(docktask.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            self.assertEqual(evaluation.can_run(), True)
            self.assertEqual(evaluation.get_error(), None)

            # evaluation task exists already
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.create_dir()
            self.assertEqual(evaluation.can_run(), False)
            self.assertEqual(evaluation.get_error(),
                             evaluation.get_dir_name() +
                             ' already exists and status is unknown')

            # evaluation task already complete
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            open(os.path.join(evaluation.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(evaluation.can_run(), False)
            self.assertEqual(evaluation.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # return immediately cause can_run is false
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.run()
            self.assertEqual(evaluation.get_error(),
                             'foo task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_evaluation_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.run()
            self.assertEqual(evaluation.get_error(),
                             'evaluation not set')
            # test files get created
            self.assertEqual(os.path.isdir(evaluation.get_dir()),
                             True)
            errfile = os.path.join(evaluation.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_pdbdb_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.evaluation = 'true'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.run()
            self.assertEqual(evaluation.get_error(),
                             'pdbdb not set')
            # test files get created
            self.assertEqual(os.path.isdir(evaluation.get_dir()),
                             True)
            errfile = os.path.join(evaluation.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_evaluation_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.evaluation = 'false'
            params.pdbdb = '/data/pdb'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.run()
            self.assertEqual(evaluation.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(evaluation.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(evaluation.get_dir(),
                                  'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(evaluation.get_dir(),
                                  'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_evaluation_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.evaluation = '/bin/doesnotexist'
            params.pdbdb = '/data/pdb'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.run()
            self.assertEqual(evaluation.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --pdbdb /data/pdb ' +
                             '--dockdir ' +
                             docktask.get_dir() + ' --outdir ' +
                             evaluation.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(evaluation.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.evaluation = 'true'
            params.pdbdb = '/data/pdb'
            docktask = D3RTask(temp_dir, params)
            docktask.set_name('foo')
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            evaluation = EvaluationTask(temp_dir, 'foo.evaluation',
                                        docktask, params)
            evaluation.run()
            self.assertEqual(evaluation.get_error(), None)
            # test files get created
            errfile = os.path.join(evaluation.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(evaluation.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(evaluation.get_dir(),
                                  'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(evaluation.get_dir(),
                                  'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
