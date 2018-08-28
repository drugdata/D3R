__author__ = 'churas'

import unittest
import tempfile
import os.path
import stat
import configparser
import httpretty

"""
test_postevaluation
--------------------------------

Tests for `postevaluation` module.
"""
import shutil
import os

from mock import Mock

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.postevaluation import PostEvaluationTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.task import SmtpEmailer
from d3r.celpp.postevaluation import PostEvaluationEmailer
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.task import WebsiteServiceConfig


class TestPostEvaluation(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = PostEvaluationTask(temp_dir, params)

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
            self.assertTrue(final_log in flist)

            # try with csv file
            csvfile = os.path.join(task.get_dir(),
                                   'Overall_RMSD_foo' +
                                   PostEvaluationTask.CSV_SUFFIX)
            open(csvfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            self.assertTrue(csvfile in flist)
            self.assertTrue(final_log in flist)

            # try with summary.txt file
            open(task.get_summary_txt(), 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            self.assertTrue(task.get_summary_txt() in flist)
            self.assertTrue(final_log in flist)

            # try with 2 csv file
            csvfile2 = os.path.join(task.get_dir(),
                                    'Overall_RMSD_foo2' +
                                    PostEvaluationTask.CSV_SUFFIX)
            open(csvfile2, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 4)
            self.assertTrue(csvfile in flist)
            self.assertTrue(csvfile2 in flist)
            self.assertTrue(task.get_summary_txt() in flist)
            self.assertTrue(final_log in flist)
        finally:
            shutil.rmtree(temp_dir)

    def test_postevaluationtask_constructor(self):
        params = D3RParameters()
        task = PostEvaluationTask('/blah', params)
        self.assertEquals(task.get_name(), PostEvaluationTask.TASK_NAME)
        self.assertEquals(task.get_stage(), 8)
        self.assertEqual(task.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(task._week_num, '0')
        self.assertEqual(task._year, '0')

        task = PostEvaluationTask('/blah/1985/dataset.week.3', params)
        self.assertEqual(task._week_num, '3')
        self.assertEqual(task._year, '1985')

        task = PostEvaluationTask(None, params)
        self.assertEqual(task._week_num, 0)
        self.assertEqual(task._year, 0)

        temp_dir = tempfile.mkdtemp()
        try:
            params.websiteserviceconfig = os.path.join(temp_dir, 'noexist')
            task = PostEvaluationTask('/blah/1985/dataset.week.3', params)
            self.assertTrue(task._webserviceconfig is not None)
        finally:
            shutil.rmtree(temp_dir)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            # post evaluation task exists already
            task = PostEvaluationTask(temp_dir, params)
            task.create_dir()
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             task.get_dir_name() +
                             ' already exists and status is unknown')

            # post evaluation task already complete
            task = PostEvaluationTask(temp_dir, params)

            open(os.path.join(task.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_get_summary_txt(self):
        params = D3RParameters()

        task = PostEvaluationTask('/foo', params)
        self.assertEqual(task.get_summary_txt(),
                         os.path.join(task.get_dir(),
                                      PostEvaluationTask.SUMMARY_TXT))

    def test_get_postevaluation_summary(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = PostEvaluationTask(temp_dir, params)
            task.create_dir()

            # no file
            res = task.get_postevaluation_summary()
            self.assertEqual(res, '\n\nNo ' + task.get_summary_txt() +
                             ' file found\n')

            # empty file
            s_file = os.path.join(task.get_dir(),
                                  PostEvaluationTask.SUMMARY_TXT)
            open(s_file, 'a').close()
            res = task.get_postevaluation_summary()
            self.assertEqual(res, '\n\n\n')

            os.chmod(task.get_summary_txt(), 0)
            res = task.get_postevaluation_summary()
            self.assertTrue('\n\nUnable to generate postevaluation' in res)

            os.chmod(task.get_summary_txt(), stat.S_IRWXU)

            # some data in the file
            f = open(s_file, 'w')
            f.write('hello\nhi\n')
            f.flush()
            f.close()
            res = task.get_postevaluation_summary()
            self.assertEqual(res, '\n\nhello\nhi\n\n')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_all_csv_files_in_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = PostEvaluationTask(temp_dir, params)
            task.create_dir()

            # test empty dir
            res = task.get_all_csv_files_in_dir()
            self.assertEqual(len(res), 0)

            # add directory with csv ending
            csvdir = os.path.join(task.get_dir(), 'hi ' +
                                  PostEvaluationTask.CSV_SUFFIX)
            os.makedirs(csvdir, mode=0o755)
            res = task.get_all_csv_files_in_dir()
            self.assertEqual(len(res), 0)

            # 1 csv
            onecsv = os.path.join(task.get_dir(), 'Overall_RMSD_hi' +
                                  PostEvaluationTask.CSV_SUFFIX)
            open(onecsv, 'a').close()
            res = task.get_all_csv_files_in_dir()
            self.assertEqual(len(res), 1)
            self.assertTrue(onecsv in res)

            # 3 csv
            twocsv = os.path.join(task.get_dir(), 'Overall_RMSD_bye' +
                                  PostEvaluationTask.CSV_SUFFIX)
            open(twocsv, 'a').close()
            threecsv = os.path.join(task.get_dir(), 'some' +
                                    PostEvaluationTask.CSV_SUFFIX)
            open(threecsv, 'a').close()

            res = task.get_all_csv_files_in_dir()
            self.assertEqual(len(res), 3)
            self.assertTrue(onecsv in res)
            self.assertTrue(twocsv in res)
            self.assertTrue(threecsv in res)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_all_evaluation_tasks(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            # no evaluation tasks
            task = PostEvaluationTask(temp_dir, params)
            res = task.get_all_evaluation_tasks()
            self.assertEqual(len(res), 0)

            # directory in path that does NOT match suffix
            task.create_dir()
            os.makedirs(os.path.join(temp_dir, 'stage.x.1'), mode=0o755)
            res = task.get_all_evaluation_tasks()
            self.assertEqual(len(res), 0)

            # one task, but is file (weird case)
            weirdtaskfile = os.path.join(temp_dir, 'stage.7.hi.' +
                                         EvaluationTaskFactory.SCORING_SUFFIX)
            open(weirdtaskfile, 'a').close()
            res = task.get_all_evaluation_tasks()
            self.assertEqual(len(res), 0)
            self.assertEqual(task.get_email_log(),
                             'Just a note, found a task with valid name, but '
                             'it is not a directory ' + weirdtaskfile)
            # one task
            etask = EvaluationTask(temp_dir, 'foo.evaluation',
                                   None, params)
            etask.create_dir()
            res = task.get_all_evaluation_tasks()
            self.assertEqual(len(res), 1)
            self.assertTrue(res[0].get_dir_name(), etask.get_dir_name())

            # three tasks two with non complete status
            etask2 = EvaluationTask(temp_dir,
                                    'bla_dock.extsubmission.evaluation',
                                    None, params)
            etask2.create_dir()
            open(os.path.join(etask2.get_dir(),
                              D3RTask.ERROR_FILE), 'a').close()

            etask3 = EvaluationTask(temp_dir,
                                    '12_x.extsubmission.evaluation',
                                    None, params)
            etask3.create_dir()
            open(os.path.join(etask3.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

            res = task.get_all_evaluation_tasks()
            self.assertEqual(len(res), 3)
            d_name_list = []
            for entry in res:
                d_name_list.append(entry.get_dir_name())

            self.assertTrue(etask.get_dir_name() in d_name_list)
            self.assertTrue(etask2.get_dir_name() in d_name_list)
            self.assertTrue(etask3.get_dir_name() in d_name_list)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_evaluationdir_args(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            # no evaluation tasks
            task = PostEvaluationTask(temp_dir, params)
            res = task._get_evaluationdir_args()
            self.assertEqual(res, '')

            # one task
            etask = EvaluationTask(temp_dir, 'foo.evaluation',
                                   None, params)
            etask.create_dir()
            res = task._get_evaluationdir_args()
            self.assertEqual(res, ' ' + PostEvaluationTask.EVALUATIONDIR_ARG +
                             ' ' + etask.get_dir())

            # three tasks
            etask2 = EvaluationTask(temp_dir,
                                    'bla_dock.extsubmission.evaluation',
                                    None, params)
            etask2.create_dir()
            open(os.path.join(etask2.get_dir(),
                              D3RTask.ERROR_FILE), 'a').close()

            etask3 = EvaluationTask(temp_dir,
                                    '12_x.extsubmission.evaluation',
                                    None, params)
            etask3.create_dir()
            cfile = os.path.join(etask3.get_dir(), D3RTask.COMPLETE_FILE)
            open(cfile, 'a').close()

            res = task._get_evaluationdir_args()
            self.assertTrue(' ' + PostEvaluationTask.EVALUATIONDIR_ARG +
                            ' ' + etask.get_dir() in res)
            self.assertTrue(' ' + PostEvaluationTask.EVALUATIONDIR_ARG +
                            ' ' + etask2.get_dir() in res)
            self.assertTrue(' ' + PostEvaluationTask.EVALUATIONDIR_ARG +
                            ' ' + etask3.get_dir() in res)
        finally:
            shutil.rmtree(temp_dir)

    def test_generate_target_object_no_blastnfilterdir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            basedir = os.path.join(temp_dir, '2000', 'dataset.week.1')
            os.makedirs(basedir, mode=0o755)
            params = D3RParameters()
            params.latest_weekly = basedir
            task = PostEvaluationTask(basedir, params)
            task.create_dir()
            res = task.generate_target_object()
            self.assertEqual(res['week'], 1)
            self.assertEqual(res['year'], 2000)
            self.assertEqual(res['portal_name'], 'notset')
            self.assertEqual(res['source'], 'notset')
            self.assertEqual(res['version'], 'unknown')
            self.assertEqual(res['targets'], 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_generate_target_object_allset(self):
        temp_dir = tempfile.mkdtemp()
        try:
            basedir = os.path.join(temp_dir, 'dataset.week')
            os.makedirs(basedir, mode=0o755)
            params = D3RParameters()
            params.latest_weekly = basedir

            btask = BlastNFilterTask(basedir, params)
            btask.create_dir()
            open(os.path.join(btask.get_dir(), '1fcy.txt'), 'a').close()
            open(os.path.join(btask.get_dir(), 'summary.txt'), 'a').close()

            cfile = os.path.join(temp_dir, 'foo.config')
            con = configparser.ConfigParser()
            con.add_section(WebsiteServiceConfig.DEFAULT)
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_PORTAL_NAME,
                    'portal')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_SOURCE, 'prod')
            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()
            params.websiteserviceconfig = cfile
            task = PostEvaluationTask(basedir, params)
            task.create_dir()
            with open(os.path.join(task.get_dir(), task.START_FILE), 'w') as f:
                f.write('4.3.2')
                f.flush()
            res = task.generate_target_object()
            self.assertEqual(res['week'], 0)
            self.assertEqual(res['year'], 0)
            self.assertEqual(res['portal_name'], 'portal')
            self.assertEqual(res['source'], 'prod')
            self.assertEqual(res['version'], '4.3.2')
            self.assertEqual(res['targets'], 1)
        finally:
            shutil.rmtree(temp_dir)

    def test_post_target_stats_to_websiteservice_targetobj_is_none(self):
        params = D3RParameters()
        task = PostEvaluationTask('/foo', params)
        task.post_target_stats_to_websiteservice(None)
        self.assertEqual(task.get_email_log(), None)

    def test_post_target_stats_to_websiteservice_config_is_none(self):
        targetobj = dict()
        targetobj['targets'] = 23
        params = D3RParameters()
        task = PostEvaluationTask('/foo', params)
        task.post_target_stats_to_websiteservice(targetobj)
        self.assertEqual(task.get_email_log(), '\nNo website service '
                                               'configuration found to '
                                               'post evaluation results\n')

    def test_post_target_stats_to_websiteservice_success_no_auth(self):
        temp_dir = tempfile.mkdtemp()
        try:
            httpretty.enable()
            httpretty.register_uri(httpretty.POST,
                                   'http://foo.com/week',
                                   body='hi')
            targetobj = dict()
            targetobj['targets'] = 23
            params = D3RParameters()
            con = configparser.ConfigParser()
            con.add_section(WebsiteServiceConfig.DEFAULT)
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_APIKEY, 'somekey')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_URL, 'http://foo.com/')
            cfile = os.path.join(temp_dir, 'web.conf')
            with open(cfile, 'w') as f:
                con.write(f)
                f.flush()
            params.websiteserviceconfig = cfile
            task = PostEvaluationTask(temp_dir, params)

            task.post_target_stats_to_websiteservice(targetobj)

            self.assertEqual(task.get_email_log(), None)
            reqy = httpretty.last_request()

            self.assertEqual(reqy.headers['X-API-KEY'], 'somekey')
            self.assertEqual(reqy.body, '{"targets": 23}')
        finally:
            shutil.rmtree(temp_dir)
            httpretty.disable()
            httpretty.reset()

    def test_post_target_stats_to_websiteservice_fail_with_auth(self):
        temp_dir = tempfile.mkdtemp()
        try:
            httpretty.enable()
            httpretty.register_uri(httpretty.POST,
                                   'http://foo.com/week',
                                   status=404,
                                   body='hi')
            targetobj = dict()
            targetobj['targets'] = 23
            params = D3RParameters()
            con = configparser.ConfigParser()
            con.add_section(WebsiteServiceConfig.DEFAULT)
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_APIKEY, 'somekey')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_URL, 'http://foo.com/')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_BASIC_USER, 'bob')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_BASIC_PASS, 'pw')

            cfile = os.path.join(temp_dir, 'web.conf')
            with open(cfile, 'w') as f:
                con.write(f)
                f.flush()
            params.websiteserviceconfig = cfile
            task = PostEvaluationTask(temp_dir, params)

            task.post_target_stats_to_websiteservice(targetobj)

            self.assertEqual(task.get_email_log(),
                             '\nwebsite REST service returned code: 404 with '
                             'text: hi and json: <bound method Response.json '
                             'of <Response [404]>>\n')
            reqy = httpretty.last_request()

            self.assertEqual(reqy.headers['X-API-KEY'], 'somekey')
            self.assertEqual(reqy.headers['Authorization'], 'Basic Ym9iOnB3')
        finally:
            shutil.rmtree(temp_dir)
            httpretty.disable()
            httpretty.reset()

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = PostEvaluationTask(temp_dir, params)
            task.create_dir()
            task.run()
            self.assertEqual(task.get_error(),
                             task.get_dir_name() +
                             ' already exists and status is ' +
                             D3RTask.UNKNOWN_STATUS)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_postevaluation_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = PostEvaluationTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(),
                             'postevaluation not set')
            # test files get created
            self.assertEqual(os.path.isdir(task.get_dir()),
                             True)
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fail_no_evaluations_and_no_emailer(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.postevaluation = 'true'
            task = PostEvaluationTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(),
                             'Did not find any Evaluation tasks to summarize')
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertTrue(os.path.isfile(os.path.join(task.get_dir(),
                                                        D3RTask.ERROR_FILE)))
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds_one_evaluation_and_no_emailer(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.postevaluation = 'echo'
            etask = EvaluationTask(temp_dir, 'foo_dock.' +
                                   EvaluationTaskFactory.SCORING_SUFFIX,
                                   None, params)
            etask.create_dir()

            task = PostEvaluationTask(temp_dir, params)
            ctask = ChallengeDataTask(temp_dir, params)
            task.run()

            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.COMPLETE_STATUS)
            cfile = os.path.join(task.get_dir(), D3RTask.COMPLETE_FILE)
            self.assertTrue(os.path.isfile(cfile))
            stdout_file = os.path.join(task.get_dir(), 'echo' +
                                       D3RTask.STDOUT_SUFFIX)
            self.assertTrue(os.path.isfile(stdout_file))

            f = open(stdout_file, 'r')
            data = f.read()
            f.close()
            ctaskdir = os.path.join(ctask.get_dir(),
                                    ctask.get_celpp_challenge_data_dir_name())

            self.assertTrue(' --evaluationdir ' + etask.get_dir() in data)
            self.assertTrue(' --challengedir ' + ctaskdir in data)

            self.assertTrue('--stageprefix stage.7. --evaluationsuffix ' +
                            '.extsubmission.evaluation$|.evaluation$' in data)

            stderr_file = os.path.join(task.get_dir(), 'echo' +
                                       D3RTask.STDERR_SUFFIX)
            self.assertTrue(os.path.isfile(stderr_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_one_evaluation_and_no_emailer(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.postevaluation = 'false'
            etask = EvaluationTask(temp_dir, 'foo_dock.' +
                                   EvaluationTaskFactory.SCORING_SUFFIX,
                                   None, params)
            etask.create_dir()

            task = PostEvaluationTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(),
                             'Non zero exit code: 1 received. '
                             'Standard out:  Standard error: ')
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertTrue(os.path.isfile(os.path.join(task.get_dir(),
                                                        D3RTask.ERROR_FILE)))
            self.assertTrue('Caught exception trying to '
                            'send evaluation email' in task.get_email_log())
            stdout_file = os.path.join(task.get_dir(), 'false' +
                                       D3RTask.STDOUT_SUFFIX)
            self.assertTrue(os.path.isfile(stdout_file))

            stderr_file = os.path.join(task.get_dir(), 'false' +
                                       D3RTask.STDERR_SUFFIX)
            self.assertTrue(os.path.isfile(stderr_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_append_to_message_log(self):
        pee = PostEvaluationEmailer('bob@bob.com', 'joe@joe.com')
        self.assertEqual(pee.get_message_log(), None)
        pee._append_to_message_log('hello')
        self.assertEqual(pee.get_message_log(), 'hello')
        pee._append_to_message_log('bye')
        self.assertEqual(pee.get_message_log(), 'hellobye')

    def test_send_postevaluation_email_none_task_and_none_email(self):
        pee = PostEvaluationEmailer('bob@bob.com', 'joe@joe.com')
        pee.send_postevaluation_email(None)
        self.assertEqual(pee.get_message_log(), '\nTask passed in is None\n')

        pee = PostEvaluationEmailer(None, None)
        pee.send_postevaluation_email('foo')
        self.assertEqual(pee.get_message_log(),
                         '\nNo email addresses in to list\n')

    def test_generate_post_evaluation_email_body(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = '1.2.3'
            weekyear = os.path.join(temp_dir, '2017', 'dataset.week.12')
            os.makedirs(weekyear, mode=0o755)
            task = PostEvaluationTask(weekyear, params)
            task.create_dir()
            f = open(task.get_summary_txt(), 'w')
            f.write('hello')
            f.flush()
            f.close()
            pee = PostEvaluationEmailer(None, None)
            subject, msg = pee._generate_post_evaluation_email_body(task)

            self.assertEqual(subject, D3RTask.SUBJECT_LINE_PREFIX +
                             '2017 week 12 post evaluation summary report')

            self.assertEqual(msg, 'Dear CELPP Admin,\n\nHere is the post '
                                  'evaluation  summary reports for CELPP '
                                  'week 12\n\n\n\nhello\n\n\nSincerely,\n\n'
                                  'CELPP Automation 1.2.3')
        finally:
            shutil.rmtree(temp_dir)

    def test_send_postevaluation_email_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            weekyear = os.path.join(temp_dir, '2017', 'dataset.week.12')
            os.makedirs(weekyear, mode=0o755)
            task = PostEvaluationTask(weekyear, params)
            task.create_dir()

            smtpemailer = SmtpEmailer()
            mockserver = D3RParameters()
            mockserver.sendmail = Mock()
            mockserver.quit = Mock()
            smtpemailer.set_alternate_smtp_server(mockserver)
            pee = PostEvaluationEmailer(['bob@bob.com'], smtpemailer)
            task.set_evaluation_emailer(pee)
            f = open(task.get_summary_txt(), 'w')
            f.write('hello')
            f.flush()
            f.close()
            pee.send_postevaluation_email(task)
            mockserver.quit.assert_any_call()
            self.assertEqual(pee.get_message_log(),
                             '\nSent post evaluation email to: bob@bob.com\n')
            self.assertEqual(mockserver.sendmail.call_count, 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_send_postevaluation_email_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            weekyear = os.path.join(temp_dir, '2017', 'dataset.week.12')
            os.makedirs(weekyear, mode=0o755)
            task = PostEvaluationTask(weekyear, params)
            task.create_dir()
            csvfile = os.path.join(task.get_dir(), 'foo.csv')
            open(csvfile, 'a').close()

            smtpemailer = SmtpEmailer()
            mockserver = D3RParameters()
            mockserver.sendmail = Mock(side_effect=IOError('ha'))
            mockserver.quit = Mock()
            smtpemailer.set_alternate_smtp_server(mockserver)
            pee = PostEvaluationEmailer(['bob@bob.com'], smtpemailer)
            task.set_evaluation_emailer(pee)
            f = open(task.get_summary_txt(), 'w')
            f.write('hello')
            f.flush()
            f.close()
            pee.send_postevaluation_email(task)
            mockserver.quit.assert_any_call()
            self.assertEqual(pee.get_message_log(),
                             '\nCaught exception trying to email :'
                             ' bob@bob.com : Caught exception ha\n')
            self.assertEqual(mockserver.sendmail.call_count, 1)
        finally:
            shutil.rmtree(temp_dir)

    # def test_run_succeeds_with_emailer(self):

    def test_run_succeeds_one_evaluation_with_emailer(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            foo_script = os.path.join(temp_dir, 'foo.py')
            f = open(foo_script, 'w')
            f.write('#! /usr/bin/env python\n')
            f.write('import sys\nimport os\n')
            f.write('outdir = sys.argv[1]\n')
            f.write('sys.stdout.write(" ".join(sys.argv))\n')
            f.write('f = open(os.path.join(outdir, "summary.txt"),"w")\n')
            f.write('f.write("summary")\n')
            f.write('f.flush()\nf.close()\n')
            f.write('f = open(os.path.join(outdir, "yo.csv"), "a")\n')
            f.write('f.write("hi")\nf.flush()\nf.close()\n')
            f.write('sys.exit(0)\n')
            f.flush()
            f.close()
            os.chmod(foo_script, stat.S_IRWXU)

            params.postevaluation = foo_script
            etask = EvaluationTask(temp_dir, 'foo_dock.' +
                                   EvaluationTaskFactory.SCORING_SUFFIX,
                                   None, params)
            etask.create_dir()
            task = PostEvaluationTask(temp_dir, params)

            # setup emailer
            smtpemailer = SmtpEmailer()
            mockserver = D3RParameters()
            mockserver.sendmail = Mock()
            mockserver.quit = Mock()
            smtpemailer.set_alternate_smtp_server(mockserver)
            pee = PostEvaluationEmailer(['bob@bob.com'], smtpemailer)
            task.set_evaluation_emailer(pee)

            ctask = ChallengeDataTask(temp_dir, params)
            task.run()

            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.COMPLETE_STATUS)
            cfile = os.path.join(task.get_dir(), D3RTask.COMPLETE_FILE)
            self.assertTrue(os.path.isfile(cfile))
            stdout_file = os.path.join(task.get_dir(), 'foo.py' +
                                       D3RTask.STDOUT_SUFFIX)
            self.assertTrue(os.path.isfile(stdout_file))

            f = open(stdout_file, 'r')
            data = f.read()
            f.close()
            ctaskdir = os.path.join(ctask.get_dir(),
                                    ctask.get_celpp_challenge_data_dir_name())

            self.assertTrue(' --evaluationdir ' + etask.get_dir() in data)
            self.assertTrue(' --challengedir ' + ctaskdir in data)

            self.assertTrue('--stageprefix stage.7. --evaluationsuffix ' +
                            '.extsubmission.evaluation$|.evaluation$' in data)

            stderr_file = os.path.join(task.get_dir(), 'foo.py' +
                                       D3RTask.STDERR_SUFFIX)

            f = open(stdout_file, 'r')
            data = f.read()
            f.close()
            self.assertTrue(' --evaluationdir ' + etask.get_dir() in data)
            self.assertTrue(os.path.isfile(stderr_file))
            self.assertEqual(task.get_error(), None)
            self.assertTrue('Sent post evaluation email to: '
                            'bob@bob.com' in task.get_email_log())
        finally:
            shutil.rmtree(temp_dir)

    """
    def test_generate_external_submission_email_body(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            dtask = D3RTask('/foo', params)
            dtask.set_name('foo')
            task = EvaluationTask(temp_dir, dtask.get_name(), dtask, params)
            task.create_dir()
            emailer = EvaluationEmailer(None, None)
            task.set_evaluation_emailer(emailer)
            rmsd = task.get_rmsd_txt()
            f = open(rmsd, 'w')
            f.write('  LMCSS  XXX\n1fcz  0.2  0.4\n')
            f.flush()
            f.close()
            subject, body = emailer\
                ._generate_external_submission_email_body(task)
            self.assertEqual(subject,
                             '[d3rcelpp] Week 0 evaluation results for foo')
            self.assertEqual(body, 'Dear CELPP Participant,\n\nHere are your '
                                   'docking evaluation results '
                                   '(RMSD, Angstroms) '
                                   'for CELPP week '
                                   '0\n\n\nEvaluation of docking'
                                   '\n=====================\n  LMCSS  XXX\n'
                                   '1fcz  0.2  0.4\n\n\n\nSincerely,\n\n'
                                   'CELPP Automation')
        finally:
            shutil.rmtree(temp_dir)

    def test_eval_emailer_append_to_message_log(self):
        emailer = EvaluationEmailer(None, None)
        self.assertEqual(emailer.get_message_log(), None)
        emailer._append_to_message_log('hi\n')
        self.assertEqual(emailer.get_message_log(), 'hi\n')
        emailer._append_to_message_log('how\n')
        self.assertEqual(emailer.get_message_log(), 'hi\nhow\n')

    def test_eval_emailer_get_external_submitter_email_none_from_guid(self):
        emailer = EvaluationEmailer(None, None)
        params = D3RParameters()
        task = EvaluationTask('/foo', None, None, params)
        self.assertEqual(emailer._get_external_submitter_email(task), None)
        self.assertEqual(emailer.get_message_log(),
                         '\nUnable to extract guid\n')


    def test_get_external_submitter_email_no_participant_email(self):

        params = D3RParameters()
        dtask = D3RTask('/foo', params)
        dtask.set_name('444' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
        task = EvaluationTask('/foo',
                              dtask.get_name(),
                              dtask, params)
        plist = [Participant('1name', '1d3rusername', '444',
                             None)]
        emailer = EvaluationEmailer(ParticipantDatabase(plist), None)
        self.assertEqual(emailer._get_external_submitter_email(task), None)
        self.assertEqual(emailer.get_message_log(),
                         '\nEmail not set for participant\n')

    def test_get_external_submitter_email_valid(self):

        params = D3RParameters()
        dtask = D3RTask('/foo', params)
        dtask.set_name('12345' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
        task = EvaluationTask('/foo',
                              dtask.get_name(),
                              dtask, params)
        plist = [Participant('1name', '1d3rusername', '12345',
                             'bob@bob.com')]
        # try single email address
        emailer = EvaluationEmailer(ParticipantDatabase(plist), None)
        emails = emailer._get_external_submitter_email(task)
        self.assertEqual(emails[0], 'bob@bob.com')
        self.assertEqual(len(emails), 1)
        self.assertEqual(emailer.get_message_log(), None)

        # try multiple email address
        plist = [Participant('1name', '1d3rusername', '12345',
                             'bob@bob.com,joe@joe.com')]
        # try single email address
        emailer = EvaluationEmailer(ParticipantDatabase(plist), None)
        emails = emailer._get_external_submitter_email(task)
        self.assertEqual(emails[0], 'bob@bob.com')
        self.assertEqual(emails[1], 'joe@joe.com')
        self.assertEqual(len(emails), 2)
        self.assertEqual(emailer.get_message_log(), None)

    def test_send_evaluation_email_none_task(self):
        emailer = EvaluationEmailer(None, None)
        emailer.send_evaluation_email(None)
        self.assertEqual(emailer.get_message_log(),
                         '\nTask passed in is None\n')

    def test_send_evaluation_email_clears_its_logs(self):
        emailer = EvaluationEmailer(None, None)
        emailer.send_evaluation_email(None)
        emailer.send_evaluation_email(None)

        self.assertEqual(emailer.get_message_log(),
                         '\nTask passed in is None\n')

    def test_send_evaluation_email_not_external_task(self):
        emailer = EvaluationEmailer(None, None)
        task = EvaluationTask('/foo', 'blah' +
                              EvaluationTask.EXT_SUBMISSION_SUFFIX,
                              None,
                              D3RParameters())
        emailer.send_evaluation_email(task)
        self.assertEqual(emailer.get_message_log(),
                         '\nNot an external submission\n')

    def test_send_evaluation_email_no_database(self):
        params = D3RParameters()
        dtask = D3RTask('/foo', params)
        dtask.set_name('444' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
        task = EvaluationTask('/foo',
                              dtask.get_name(),
                              dtask, params)
        emailer = EvaluationEmailer(None, None)
        emailer.send_evaluation_email(task)
        self.assertEqual(emailer.get_message_log(),
                         '\nParticipant database is None cannot send docking '
                         'evaluation email!!!\n')

    def test_send_external_submission_email_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.program = 'foo'
            params.version = '1'
            dtask = D3RTask('/foo', params)
            dtask.set_name('12345' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
            task = EvaluationTask(temp_dir,
                                  dtask.get_name(),
                                  dtask, params)
            task.create_dir()
            f = open(task.get_rmsd_txt(), 'w')
            f.write('hi\n')
            f.flush()
            f.close()
            plist = [Participant('1name', '1d3rusername', '12345',
                                 'bob@bob.com')]
            # try single email address
            smtpemailer = SmtpEmailer()
            mockserver = D3RParameters()
            mockserver.sendmail = Mock()
            mockserver.quit = Mock()
            smtpemailer.set_alternate_smtp_server(mockserver)
            emailer = EvaluationEmailer(ParticipantDatabase(plist), None)

            emailer.set_alternate_smtp_emailer(smtpemailer)
            emailer.send_evaluation_email(task)
            mockserver.quit.assert_any_call()
            self.assertEqual(emailer.get_message_log(),
                             '\nSent evaluation email to: bob@bob.com\n')
            self.assertEqual(mockserver.sendmail.call_count, 1)

        finally:
            shutil.rmtree(temp_dir)

    def test_send_external_submission_email_sendmail_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.program = 'foo'
            params.version = '1'
            dtask = D3RTask('/foo', params)
            dtask.set_name('12345' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
            task = EvaluationTask(temp_dir,
                                  dtask.get_name(),
                                  dtask, params)
            plist = [Participant('1name', '1d3rusername', '12345',
                                 'bob@bob.com')]
            # try single email address
            smtpemailer = SmtpEmailer()
            mockserver = D3RParameters()
            mockserver.sendmail = Mock(side_effect=IOError('ha'))
            mockserver.quit = Mock()
            smtpemailer.set_alternate_smtp_server(mockserver)
            emailer = EvaluationEmailer(ParticipantDatabase(plist), None)
            emailer.set_alternate_smtp_emailer(smtpemailer)
            emailer.send_evaluation_email(task)
            mockserver.quit.assert_any_call()
            self.assertEqual(emailer.get_message_log(),
                             '\nCaught exception trying to email '
                             'participant : Caught exception ha\n')

        finally:
            shutil.rmtree(temp_dir)

    def test_send_external_submission_email_no_submitter_email(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.program = 'foo'
            params.version = '1'
            dtask = D3RTask('/foo', params)
            dtask.set_name('12345' + EvaluationTask.EXT_SUBMISSION_SUFFIX)
            task = EvaluationTask(temp_dir,
                                  dtask.get_name(),
                                  dtask, params)
            plist = [Participant('1name', '1d3rusername', '1234',
                                 'bob@bob.com')]
            # try single email address
            emailer = EvaluationEmailer(ParticipantDatabase(plist), None)
            smtpemailer = SmtpEmailer()
            emailer.set_alternate_smtp_emailer(smtpemailer)
            emailer.send_evaluation_email(task)
            self.assertEqual(emailer.get_message_log(),
                             '\nNo participant found with guid: 12345\n')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_reply_to_address(self):
        # test get reply_to where no replytoaddress is set
        emailer = EvaluationEmailer(None, None)
        val = emailer._get_reply_to_address('bob@bob.com')
        self.assertEqual(val, 'bob@bob.com')

        emailer = EvaluationEmailer(None, 'joe@joe.com')
        val = emailer._get_reply_to_address('bob@bob.com')
        self.assertEqual(val, 'joe@joe.com')

    def test_get_smtp_emailer_valid(self):
        try:
            emailer = EvaluationEmailer(None, None)
            ss = emailer._get_smtp_emailer()
            self.assertTrue(ss is not None)
        finally:
            pass
    """

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
