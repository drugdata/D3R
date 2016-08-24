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
from d3r.celpp.dataimport import DataImportTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.glide import GlideTask
from d3r.celpp.evaluation import PathNotDirectoryError


class TestEvaluation(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = EvaluationTask(temp_dir, 'glide',
                                  GlideTask(temp_dir, params), params)
            # try with no dir
            self.assertEqual(task.get_uploadable_files(), [])

            # try with empty dir
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # try with final log
            final_log = os.path.join(task.get_dir(),
                                     EvaluationTask.FINAL_LOG)
            open(final_log, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(final_log)

            # try with RMSD.txt
            rmsd = os.path.join(task.get_dir(),
                                EvaluationTask.RMSD_TXT)
            open(rmsd, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(rmsd)

            # try with empty pbdid dir
            pbdid = os.path.join(task.get_dir(), '8www')
            os.mkdir(pbdid)
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(rmsd)

            # try with score/rot-LMCSS_doc_pv_complex1.pdb
            score = os.path.join(pbdid, 'score')
            os.mkdir(score)
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)

            LMCSS = os.path.join(score, 'rot-LMCSS_dock_pv_complex1.pdb')
            open(LMCSS, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(LMCSS)

            # try with score/rot-SMCSS_doc_pv_complex1.pdb
            SMCSS = os.path.join(score, 'rot-SMCSS_dock_pv_complex1.pdb')
            open(SMCSS, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 4)
            flist.index(SMCSS)

            # try with score/rot-hiResApo_doc_pv_complex1.pdb
            hiResApo = os.path.join(score, 'rot-hiResApo_dock_pv_complex1.pdb')
            open(hiResApo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 5)
            flist.index(hiResApo)

            # try with score/rot-hiResHolo_doc_pv_complex1.pdb
            hiResHolo = os.path.join(score,
                                     'rot-hiResHolo_dock_pv_complex1.pdb')
            open(hiResHolo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 6)
            flist.index(hiResHolo)

            # try with score/crystal.pdb
            crystal = os.path.join(score, 'crystal.pdb')
            open(crystal, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 7)
            flist.index(crystal)

            # try with stderr/stdout files
            errfile = os.path.join(task.get_dir(), 'evaluate.py.stderr')
            open(errfile, 'a').close()
            outfile = os.path.join(task.get_dir(), 'evaluate.py.stdout')
            open(outfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 9)
            flist.index(crystal)
            flist.index(hiResHolo)
            flist.index(hiResApo)
            flist.index(SMCSS)
            flist.index(LMCSS)
            flist.index(errfile)
            flist.index(outfile)
            flist.index(final_log)
            flist.index(rmsd)
        finally:
            shutil.rmtree(temp_dir)

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
            dataimport = DataImportTask(temp_dir, params)
            blast = BlastNFilterTask(temp_dir, params)
            os.mkdir(os.path.join(temp_dir, dataimport.get_dir_name()))
            os.mkdir(os.path.join(temp_dir, blast.get_dir_name()))

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
                                  EvaluationTaskFactory.DOCKSTAGE_PREFIX +
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
                                    EvaluationTaskFactory.DOCKSTAGE_PREFIX +
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
                                    EvaluationTaskFactory.DOCKSTAGE_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.COMPLETE_FILE), 'a').close()

            freddir = os.path.join(temp_dir,
                                   EvaluationTaskFactory.DOCKSTAGE_PREFIX +
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
                                    EvaluationTaskFactory.DOCKSTAGE_PREFIX +
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
                                    EvaluationTaskFactory.DOCKSTAGE_PREFIX +
                                    'glide')
            os.mkdir(glidedir)
            open(os.path.join(glidedir, D3RTask.ERROR_FILE), 'a').close()
            freddir = os.path.join(temp_dir,
                                   EvaluationTaskFactory.DOCKSTAGE_PREFIX +
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
        self.assertEquals(evaluation.get_stage(), 7)

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

    def test_get_evaluation_summary(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = EvaluationTask(temp_dir, 'foo', None, params)
            # test where no RMSD.txt file is found
            val = task._get_evaluation_summary()
            self.assertEqual(val,
                             '\nEvaluation of docking\n'
                             '=====================\nNo ' +
                             task.get_rmsd_txt() + ' file found.\n')

            # test with valid RMSD.txt file
            task.create_dir()
            f = open(task.get_rmsd_txt(), 'w')
            f.write('LMCSS\n1fcz 0.465\n')
            f.flush()
            f.close()
            val = task._get_evaluation_summary()
            self.assertEqual(val,
                             '\nEvaluation of docking\n'
                             '=====================\n'
                             'LMCSS\n1fcz 0.465\n\n')

            # test where reading RMSD.txt throws exception
            os.chmod(task.get_rmsd_txt(), 0)
            val = task._get_evaluation_summary()
            self.assertTrue(val.startswith('\nEvaluation of docking\n'
                                           '=====================\n'
                                           'Unable to generate evaluation'
                                           ' summary ('))
        finally:
            shutil.rmtree(temp_dir)

    def test_get_rmsd_txt(self):
        params = D3RParameters()
        task = EvaluationTask('/ha', 'foo', None, params)
        self.assertEqual(task.get_rmsd_txt(),
                         os.path.join('/ha', task.get_dir_name(),
                                      EvaluationTask.RMSD_TXT))

    def test_am_i_an_external_submission(self):
        params = D3RParameters()
        task = EvaluationTask('/ha', 'foo', None, params)
        self.assertEqual(task._am_i_an_external_submission(), False)
        task.set_name('blah' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
        self.assertEqual(task._am_i_an_external_submission(), True)

    def test_get_participant_database_no_csv_file_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            etf = EvaluationTaskFactory(temp_dir, params)
            pdb = etf._get_participant_database()
            self.assertEqual(pdb, None)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_participant_database_with_valid_csvfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            dimport = DataImportTask(temp_dir, params)
            dimport.create_dir()
            csvfile = dimport.get_participant_list_csv()
            f = open(csvfile, 'w')
            f.write('name,d3rusername,guid,email\n')
            f.write('joe,jj,123,j@j.com\n')
            f.flush()
            f.close()
            etf = EvaluationTaskFactory(temp_dir, params)
            pdb = etf._get_participant_database()
            p = pdb.get_participant_by_guid('123')
            self.assertEqual(p.get_email(), 'j@j.com')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_guid_from_task_name_for_external_task(self):
        params = D3RParameters()
        task = EvaluationTask('/foo', '123' +
                              EvaluationTask.EXT_SUBMISSION_SUFFIX,
                              None,
                              params)
        self.assertEqual(task._get_guid_from_task_name_for_external_task(),
                         '123')

    def test_send_external_submission_email(self):
        self.assertEqual(1, 2, 'Need to finish implementation and write tests')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
