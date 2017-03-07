#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_celpprunner
----------------------------------

Tests for `celpprunner` module.
"""

import unittest
import tempfile
import os
import os.path
import shutil
import gzip
from datetime import date

from d3r import celpprunner
from d3r.celpp.task import D3RParameters
from d3r.celpp import util
from d3r.celpp.task import D3RTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp.proteinligprep import ProteinLigPrepTask
from d3r.celpp.glide import GlideTask
from d3r.celpp.vina import AutoDockVinaTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.chimeraprep import ChimeraProteinLigPrepTask


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
    """Tests celpprunner command line script
    """
    param = D3RParameters()

    blast = BlastNFilterTask('/foo', param)
    BLAST_DIR_NAME = blast.get_dir_name()
    BLAST_NAME = blast.get_name()

    data = DataImportTask('/foo', param)
    IMPORT_DIR_NAME = data.get_dir_name()
    IMPORT_NAME = data.get_name()

    makedb = MakeBlastDBTask('/foo', param)
    MAKEDB_DIR_NAME = makedb.get_dir_name()
    MAKEDB_NAME = makedb.get_name()

    glide = GlideTask('/foo', param)
    GLIDE_DIR_NAME = glide.get_dir_name()

    prot = ProteinLigPrepTask('/foo', param)
    PROT_DIR_NAME = prot.get_dir_name()

    vina = AutoDockVinaTask('/foo', param)
    VINA_DIR_NAME = vina.get_dir_name()

    chall = ChallengeDataTask('/foo', param)
    CHALL_DIR_NAME = chall.get_dir_name()
    CHALL_NAME = chall.get_name()

    chimeraprep = ChimeraProteinLigPrepTask('/foo', param)
    CHIMERAPREP_DIR_NAME = chimeraprep.get_dir_name()

    def setUp(self):
        pass

    def test_get_lock(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.latest_weekly = temp_dir

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

    def test_parse_arguments(self):
        theargs = ['--stage', 'blast', 'foo']
        result = celpprunner._parse_arguments('hi', theargs)
        self.assertEqual(result.stage, 'blast')
        self.assertEqual(result.celppdir, 'foo')
        self.assertEqual(result.email, None)
        self.assertEqual(result.summaryemail, None)
        self.assertEqual(result.loglevel, celpprunner.DEFAULT_LOG_LEVEL)
        self.assertEqual(result.blastnfilter, 'blastnfilter.py')
        self.assertEqual(result.proteinligprep, 'proteinligprep.py')
        self.assertEqual(result.evaluation, 'evaluate.py')
        self.assertEqual(result.makeblastdb, 'makeblastdb')
        self.assertEqual(result.genchallenge, 'genchallengedata.py')
        self.assertEqual(result.chimeraprep, 'chimera_proteinligprep.py')
        self.assertEqual(result.skipimportwait, False)
        self.assertEqual(result.importretry, 60)
        self.assertEqual(result.importsleep, 600)
        self.assertEqual(result.rdkitpython, '')
        theargs = ['foo', '--stage', 'dock,glide', '--email', 'b@b.com,h@h',
                   '--log', 'ERROR',
                   '--blastnfilter', '/bin/blastnfilter.py',
                   '--proteinligprep', '/bin/proteinligprep.py',
                   '--postanalysis', '/bin/postanalysis.py',
                   '--glide', '/bin/glide.py',
                   '--vina', '/bin/vina.py',
                   '--customweekdir',
                   '--evaluation', '/bin/evaluation.py',
                   '--makeblastdb', '/bin/makeblastdb',
                   '--genchallenge', '/bin/gen.py',
                   '--chimeraprep', '/bin/chimeraprep.py',
                   '--skipimportwait',
                   '--importretry', '10',
                   '--importsleep', '30',
                   '--rdkitpython', '/usr/bin',
                   '--summaryemail', 'j@j,g@g']
        result = celpprunner._parse_arguments('hi', theargs)
        self.assertEqual(result.stage, 'dock,glide')
        self.assertEqual(result.celppdir, 'foo')
        self.assertEqual(result.email, 'b@b.com,h@h')
        self.assertEqual(result.summaryemail, 'j@j,g@g')
        self.assertEqual(result.loglevel, 'ERROR')
        self.assertEqual(result.blastnfilter, '/bin/blastnfilter.py')
        self.assertEqual(result.proteinligprep, '/bin/proteinligprep.py')
        self.assertEquals(result.postanalysis, '/bin/postanalysis.py')
        self.assertEquals(result.glide, '/bin/glide.py')
        self.assertEquals(result.evaluation, '/bin/evaluation.py')
        self.assertEquals(result.customweekdir, True)
        self.assertEqual(result.makeblastdb, '/bin/makeblastdb')
        self.assertEqual(result.vina, '/bin/vina.py')
        self.assertEqual(result.genchallenge, '/bin/gen.py')
        self.assertEqual(result.chimeraprep, '/bin/chimeraprep.py')
        self.assertEqual(result.skipimportwait, True)
        self.assertEqual(result.importretry, 10)
        self.assertEqual(result.importsleep, 30)
        self.assertEqual(result.rdkitpython, '/usr/bin')

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
        self.assertEquals(task_two._run_count, 1)

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
                          os.path.join('foo', TestCelppRunner.BLAST_DIR_NAME))

        task_list = celpprunner.get_task_list_for_stage(params,
                                                        'proteinligprep')
        self.assertEquals(len(task_list), 1)

        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', TestCelppRunner.PROT_DIR_NAME))

        task_list = celpprunner.get_task_list_for_stage(params, 'import')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', TestCelppRunner.IMPORT_DIR_NAME))

        task_list = celpprunner.get_task_list_for_stage(params, 'glide')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', TestCelppRunner.GLIDE_DIR_NAME))

        task_list = celpprunner.get_task_list_for_stage(params, 'vina')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', TestCelppRunner.VINA_DIR_NAME))

        task_list = celpprunner.get_task_list_for_stage(params,
                                                        'challengedata')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo', TestCelppRunner.CHALL_DIR_NAME))

        task_list = celpprunner.get_task_list_for_stage(params,
                                                        'chimeraprep')
        self.assertEquals(len(task_list), 1)
        self.assertEquals(task_list[0].get_dir(),
                          os.path.join('foo',
                                       TestCelppRunner.CHIMERAPREP_DIR_NAME))

    def test_get_task_list_for_stage_createchallenge(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.latest_weekly = temp_dir
            task_list = celpprunner.get_task_list_for_stage(
                params, celpprunner.CREATE_CHALLENGE)
            self.assertEqual(len(task_list), 4)
            self.assertEqual(task_list[0].get_name(),
                             TestCelppRunner.MAKEDB_NAME)
            self.assertEqual(task_list[1].get_name(),
                             TestCelppRunner.IMPORT_NAME)
            self.assertEqual(task_list[2].get_name(),
                             TestCelppRunner.BLAST_NAME)
            self.assertEqual(task_list[3].get_name(),
                             TestCelppRunner.CHALL_NAME)

        finally:
            shutil.rmtree(temp_dir)

    def test_get_task_list_for_stage_for_scoring_stage_with_nonefound(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.latest_weekly = temp_dir
            try:
                celpprunner.get_task_list_for_stage(params, 'evaluation')
            except NotImplementedError as e:
                self.assertEqual(e.message,
                                 'uh oh no tasks for evaluation stage')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_task_list_for_stage_for_scoring_stage_with_onefound(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.latest_weekly = temp_dir
            glidedir = os.path.join(temp_dir,
                                    EvaluationTaskFactory.DOCKSTAGE_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.COMPLETE_FILE), 'a').close()
            task_list = celpprunner.get_task_list_for_stage(params,
                                                            'evaluation')
            self.assertEqual(len(task_list), 1)
            self.assertEqual(task_list[0].get_name(), 'glide.evaluation')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_task_list_for_stage_for_scoring_stage_with_twofound(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.latest_weekly = temp_dir
            glidedir = os.path.join(temp_dir,
                                    EvaluationTaskFactory.DOCKSTAGE_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.COMPLETE_FILE), 'a').close()
            freddir = os.path.join(temp_dir,
                                   EvaluationTaskFactory.DOCKSTAGE_PREFIX +
                                   'fred')
            os.mkdir(freddir)
            open(os.path.join(freddir, D3RTask.COMPLETE_FILE), 'a').close()

            task_list = celpprunner.get_task_list_for_stage(params,
                                                            'evaluation')
            self.assertEqual(len(task_list), 2)
        finally:
            shutil.rmtree(temp_dir)

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
            makedb_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                      TestCelppRunner.MAKEDB_DIR_NAME)
            os.makedirs(makedb_dir)
            open(os.path.join(makedb_dir, 'complete'), 'a').close()
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

            makedb_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                      TestCelppRunner.MAKEDB_DIR_NAME)
            os.makedirs(makedb_dir)
            open(os.path.join(makedb_dir, 'complete'), 'a').close()

            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        TestCelppRunner.IMPORT_DIR_NAME)
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
            os.mkdir(os.path.join(temp_dir, TestCelppRunner.MAKEDB_DIR_NAME))
            open(os.path.join(temp_dir, TestCelppRunner.MAKEDB_DIR_NAME,
                              'error'),
                 'a').close()
            os.mkdir(os.path.join(temp_dir, '2015'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.1'))
            self.assertEqual(celpprunner.run_stages(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_challenge_and_proteinligprep_no_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.pdbdb = '/pdbdb'
            theargs.celppdir = os.path.join(temp_dir)
            theargs.stage = 'challengedata,proteinligprep'

            blastdb_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                       TestCelppRunner.BLAST_DIR_NAME)
            os.makedirs(blastdb_dir)
            open(os.path.join(blastdb_dir, 'complete'), 'a').close()

            theargs.proteinligprep = 'echo'
            theargs.genchallenge = 'echo'

            self.assertEqual(celpprunner.run_stages(theargs), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_blast_and_proteinligprep_blast_has_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir)
            theargs.stage = 'blast,proteinligprep'
            os.mkdir(os.path.join(temp_dir, TestCelppRunner.MAKEDB_DIR_NAME))
            open(os.path.join(temp_dir, TestCelppRunner.MAKEDB_DIR_NAME,
                              'complete'),
                 'a').close()
            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        TestCelppRunner.IMPORT_DIR_NAME)
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir, D3RTask.ERROR_FILE), 'a').close()
            theargs.blastnfilter = 'echo'
            theargs.proteinligprep = 'echo'
            self.assertEqual(celpprunner.run_stages(theargs), 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_makedb_blast_chall_proteinligprep_glide_no_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.pdbdb = '/pdbdb'
            theargs.celppdir = os.path.join(temp_dir)

            theargs.stage = 'makedb,blast,challengedata,proteinligprep,glide'

            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        TestCelppRunner.IMPORT_DIR_NAME)
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir,
                              D3RTask.COMPLETE_FILE), 'a').close()

            fakegz = os.path.join(temp_dir, 'fake.gz')

            f = gzip.open(fakegz, 'wb')
            f.write('hello\n')
            f.flush()
            f.close()

            theargs.pdbsequrl = 'file://'+fakegz

            theargs.makeblastdb = 'echo'
            theargs.blastnfilter = 'echo'
            theargs.postanalysis = 'true'
            theargs.proteinligprep = 'echo'
            theargs.glide = 'echo'
            theargs.genchallenge = 'echo'
            self.assertEqual(celpprunner.run_stages(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_stages_makedb_through_glide(self):
        """This should test the following stages will run
           makedb,import,blast,challengedata,proteinligprep,glide,vina
        """
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = D3RParameters()
            theargs.pdbdb = '/pdbdb'
            theargs.celppdir = os.path.join(temp_dir)

            theargs.stage = 'makedb,import,blast,challengedata,proteinligprep,' \
                            'chimeraprep,glide,vina'

            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        TestCelppRunner.IMPORT_DIR_NAME)
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir,
                              D3RTask.COMPLETE_FILE), 'a').close()

            fakegz = os.path.join(temp_dir, 'fake.gz')

            f = gzip.open(fakegz, 'wb')
            f.write('hello\n')
            f.flush()
            f.close()

            theargs.pdbsequrl = 'file://' + fakegz
            theargs.pdbfileurl = 'file://' + fakegz

            theargs.compinchi = 'file://' + fakegz
            theargs.version = '1.0.0'
            theargs.makeblastdb = 'echo'
            theargs.blastnfilter = 'echo'
            theargs.postanalysis = 'true'
            theargs.proteinligprep = 'echo'
            theargs.glide = 'echo'
            theargs.vina = 'echo'
            theargs.genchallenge = 'echo'
            theargs.chimeraprep = 'echo'
            self.assertEqual(celpprunner.run_stages(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_get_task_list_for_stage_extsubmission(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.pdbdb = '/pdbdb'
            theargs.latest_weekly = temp_dir
            theargs.stage = 'extsubmission'
            try:
                celpprunner.get_task_list_for_stage(theargs, 'extsubmission')
                self.fail('expected NotImplementedError')
            except NotImplementedError:
                pass

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

    def test_run_stages_customweekdir_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = temp_dir
            theargs.customweekdir = True
            theargs.createweekdir = True
            theargs.stage = ''
            try:
                self.assertEquals(celpprunner.run_stages(theargs), 0)
                self.fail('Expected NotImplementedError')
            except NotImplementedError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_main_success(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = ['celpprunner.py', '--stage',
                       'blast', '--pdbdb', '/pdbdb',
                       '--blastnfilter', 'echo',
                       '--postanalysis', 'true',
                       temp_dir]

            makedb_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                      TestCelppRunner.MAKEDB_DIR_NAME)
            os.makedirs(makedb_dir)
            open(os.path.join(makedb_dir, 'complete'), 'a').close()

            d_import_dir = os.path.join(temp_dir, '2015', 'dataset.week.1',
                                        TestCelppRunner.IMPORT_DIR_NAME)
            os.makedirs(d_import_dir)
            open(os.path.join(d_import_dir, 'complete'), 'a').close()

            self.assertEqual(celpprunner.main(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_main_where_run_stages_raises_error(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = ['celpprunner.py', '--stage',
                       'foo', os.path.join(temp_dir, 'notexistdir')]
            self.assertEqual(celpprunner.main(theargs), 2)

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
