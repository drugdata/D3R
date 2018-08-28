__author__ = 'churas'

import os
import logging
import requests
from requests.auth import HTTPBasicAuth
from d3r.celpp.task import Attachment
from d3r.celpp import util
from d3r.celpp.task import D3RTask
from d3r.celpp.task import WebsiteServiceConfig
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.blastnfilter import BlastNFilterTask

logger = logging.getLogger(__name__)


class PathNotDirectoryError(Exception):
    """Path is not a directory Error
    """
    pass


class PostEvaluationEmailer(object):
    """Sends evaluation email
    """
    def __init__(self, to_list, emailer):
        self._emailer = emailer
        self._msg_log = None
        self._to_list = to_list

    def _append_to_message_log(self, msg):
        if self._msg_log is None:
            self._msg_log = str(msg)
            return
        self._msg_log += str(msg)

    def get_message_log(self):
        """Gets the error log for this object
        :returns: Error log as string
        """
        return self._msg_log

    def _generate_post_evaluation_email_body(self, petask):
        """Creates body of email and subject
        :returns subject,body: as strings
        """
        weekno = util.get_celpp_week_number_from_path(petask.get_path())
        year = util.get_celpp_year_from_path(petask.get_path())

        msg = ('Dear CELPP Admin,\n\nHere is the post evaluation  ' +
               'summary reports for CELPP week ' + str(weekno) + '\n\n')

        msg += petask.get_postevaluation_summary()

        msg += ('\n\nSincerely,\n\nCELPP Automation ' +
                petask.get_program_version())

        subject = (D3RTask.SUBJECT_LINE_PREFIX + str(year) + ' ' +
                   'week ' + str(weekno) + ' post evaluation summary report')
        return subject, msg

    def send_postevaluation_email(self, petask):
        """Sends evaluation email appending issues to message log
        """

        # clear message log
        self._msg_log = None
        if petask is None:
            logger.error('Task passed in is None')
            self._append_to_message_log('\nTask passed in is None\n')
            return

        try:

            if self._to_list is None:
                logger.debug('No email addresses in to list')
                self._append_to_message_log('\nNo email addresses ' +
                                            'in to list\n')
                return

            subject, msg = self._generate_post_evaluation_email_body(petask)

            attach_list = None
            for csv_file in petask.get_all_csv_files_in_dir():
                if attach_list is None:
                    attach_list = []
                attach_list.append(Attachment(csv_file,
                                              os.path.basename(csv_file)))

            self._emailer.send_email(self._to_list, subject, msg,
                                     attachments=attach_list)

            self._append_to_message_log('\nSent post evaluation email to: ' +
                                        ", ".join(self._to_list) + '\n')
        except Exception as e:
            logger.exception('Caught exception')
            self._append_to_message_log('\nCaught exception trying to email '
                                        ': ' +
                                        ', '.join(self._to_list) + ' : ' +
                                        str(e) + '\n')


class PostEvaluationTask(D3RTask):
    """Generates summary of all evaluations run
    """
    RMSD_TXT = 'RMSD.txt'
    RMSD_PICKLE = 'RMSD.pickle'
    TASK_NAME = 'postevaluation'
    FINAL_LOG = 'final.log'
    CSV_SUFFIX = '.csv'
    SUMMARY_TXT = 'summary.txt'
    EVALUATIONDIR_ARG = '--evaluationdir'

    def __init__(self, path, args):
        super(PostEvaluationTask, self).__init__(path, args)
        self.set_name(PostEvaluationTask.TASK_NAME)
        et = EvaluationTask('/foo', '', None, None)
        self._eval_task_prefix_str = et.get_dir_name()
        self._eval_task_suffix_str = ('' +
                                      EvaluationTask.EXT_SUBMISSION_SUFFIX +
                                      '.' +
                                      EvaluationTaskFactory.SCORING_SUFFIX +
                                      '$|.' +
                                      EvaluationTaskFactory.SCORING_SUFFIX +
                                      '$')
        self.set_stage(et.get_stage() + 1)
        self._emailer = None
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._webserviceconfig = None

        try:
            config = self.get_args().websiteserviceconfig
            self._webserviceconfig = WebsiteServiceConfig(configfile=config)
        except AttributeError:
            pass
        try:
            self._week_num = util.\
                get_celpp_week_number_from_path(self.get_path())
            self._year = util.get_celpp_year_from_path(self.get_dir())
        except Exception:
            self._week_num = 0
            self._year = 0

    def set_evaluation_emailer(self, emailer):
        self._emailer = emailer

    def get_all_evaluation_tasks(self):
        """Look for all  evaluation tasks
        """
        task_list = []
        prefix_len = len(self._eval_task_prefix_str)
        base_dir = self.get_path()
        for entry in os.listdir(base_dir):
            if not entry.endswith('.' + EvaluationTaskFactory.SCORING_SUFFIX):
                continue
            full_path = os.path.join(base_dir, entry)
            if not os.path.isdir(full_path):
                logger.info('Suffix matched, but is not directory: ' +
                            full_path)
                self.append_to_email_log('Just a note, found a task ' +
                                         'with valid name, but it is not '
                                         'a directory ' + full_path)
                continue
            task = EvaluationTask(base_dir, entry[prefix_len:], None,
                                  self.get_args())
            task.update_status_from_filesystem()
            if task.get_status() != D3RTask.COMPLETE_STATUS:
                logger.info('fyi ' + task.get_name() + ' task ' +
                            ' has a status of ' + task.get_status())
            task_list.append(task)
        return task_list

    def _get_evaluationdir_args(self):
        """Generates a string containing all the evaluation
           tasks to pass to post_evaluation.py script
        :returns: string of format --evaluationdir <path1>
                  --evaluationdir <path2>
        """
        eval_args = ''
        for task in self.get_all_evaluation_tasks():
            logger.debug('Adding ' + task.get_name() + ' to arglist')
            eval_args += (' ' + PostEvaluationTask.EVALUATIONDIR_ARG + ' ' +
                          task.get_dir())
        return eval_args

    def get_all_csv_files_in_dir(self):
        """Gets all the files ending with CSV_SUFFIX
           in task directory
        """
        csv_list = []
        out_dir = self.get_dir()
        for entry in os.listdir(out_dir):
            if entry.endswith(PostEvaluationTask.CSV_SUFFIX):
                full_path = os.path.join(out_dir, entry)
                if os.path.isfile(full_path):
                    csv_list.append(full_path)
                else:
                    logger.warning(full_path +
                                   ' is a directory which is weird')
        return csv_list

    def get_summary_txt(self):
        """Gets summary.txt file path
        """
        return os.path.join(self.get_dir(), PostEvaluationTask.SUMMARY_TXT)

    def get_postevaluation_summary(self):
        """Summarizes post evaluation results into a string

        """
        summary_file = self.get_summary_txt()
        start_line = '\n\n'
        if not os.path.isfile(summary_file):
            return start_line + 'No ' + summary_file + ' file found\n'

        try:
            f = open(summary_file, 'r')
            summary = f.read()
            f.close()
            return start_line + summary + '\n'
        except Exception as e:
            logger.exception('Caught exception')
            return (start_line + 'Unable to generate postevaluation '
                                 'summary (' + str(e) + ')\n')

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files if found on the file system
           plus stderr/stdout files

           *.csv
           final.log
           :returns: list of files that can be uploaded
        """
        # get the stderr/stdout files
        file_list = super(PostEvaluationTask, self).get_uploadable_files()
        out_dir = self.get_dir()
        try:
            final_log = os.path.join(out_dir, PostEvaluationTask.FINAL_LOG)
            if os.path.isfile(final_log):
                file_list.append(final_log)

            if os.path.isfile(self.get_summary_txt()):
                file_list.append(self.get_summary_txt())

            csv_files = self.get_all_csv_files_in_dir()
            logger.debug('Appending ' + str(len(csv_files)) + ' csv files')
            file_list.extend(csv_files)
        except OSError:
            logger.exception('Caught exception looking for csv files')
        return file_list

    def can_run(self):
        """Determines if task can actually run

           The method verifies this task does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        # check this task is not complete and does not exist

        self.update_status_from_filesystem()
        if self.get_status() == D3RTask.COMPLETE_STATUS:
            logger.debug("No work needed for " + self.get_name() +
                         " task")
            return False

        if self.get_status() != D3RTask.NOTFOUND_STATUS:
            logger.warning(self.get_name() + " task was already " +
                           "attempted, but there was a problem")
            self.set_error(self.get_dir_name() + ' already exists and ' +
                           'status is ' + self.get_status())
            return False
        self._can_run = True
        return True

    def generate_target_object(self):
        blast = BlastNFilterTask(self.get_args().latest_weekly,
                                 self.get_args())

        summary = blast.get_blastnfilter_summary()

        numtargets = summary.get_targets_found()
        targetobj = {}
        version = self.get_version_from_start_file()

        wsource = 'notset'
        pname = 'notset'
        if self._webserviceconfig is not None:
            wsource = self._webserviceconfig.get_source()
            pname = self._webserviceconfig.get_portal_name()

        targetobj['version'] = version
        targetobj['source'] = wsource
        targetobj['week'] = int(self._week_num)
        targetobj['year'] = int(self._year)
        targetobj['portal_name'] = pname
        targetobj['targets'] = numtargets
        return targetobj

    def post_target_stats_to_websiteservice(self, targetobj):
        """Posts target statistices to website REST service
        """
        if targetobj is None:
            logger.debug('targetobj is None. Skipping post to service')
            return

        if self._webserviceconfig is None:
            logger.warning('No service information available to post '
                           'evaluation results')
            self.append_to_email_log('\nNo website service configuration '
                                     'found to post evaluation results\n')
            return

        theheader = dict()

        theheader[WebsiteServiceConfig.API_KEY_KEY] = self.\
            _webserviceconfig.get_apikey()
        theheader[WebsiteServiceConfig.
                  CONTENT_TYPE_KEY] = WebsiteServiceConfig.CONTENT_TYPE_VAL
        bauth = None
        if self._webserviceconfig.get_basicauth_user() is not None:
            bauth = HTTPBasicAuth(self._webserviceconfig.get_basicauth_user(),
                                  self._webserviceconfig.
                                  get_basicauth_password())
        logger.debug('Posting targetobj to ' +
                     str(self._webserviceconfig.get_targets_url()))
        r = requests.post(self._webserviceconfig.get_targets_url(),
                          headers=theheader, json=targetobj, auth=bauth,
                          timeout=self._webserviceconfig.get_timeout())

        if r.status_code != 200:
            msg = ('website REST service returned code: ' +
                   str(r.status_code) + ' with text: ' + str(r.text) +
                   ' and json: ' + str(r.json))
            logger.error(msg)
            self.append_to_email_log('\n' + msg + '\n')
        else:
            logger.debug('Post returned code 200. Success')
        return

    def run(self):
        """Runs PostEvaluationTask

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes evaluation script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(PostEvaluationTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('postevaluation set to ' +
                         self.get_args().postevaluation)
        except AttributeError:
            self.set_error('postevaluation not set')
            self.end()
            return

        #
        # --evaluationdir <path to stage.X.evaluation>
        # --evaluationdir <path to another etc...>
        # --challengedir <path to challenge dir ie
        #                 stage.4.challengedata/celpp_week10_2017>
        #
        chall = ChallengeDataTask(self.get_path(), self.get_args())
        chall_dname = chall.get_celpp_challenge_data_dir_name()

        evaldir_args = self._get_evaluationdir_args()
        if evaldir_args is not '':
            cmd_to_run = (self.get_args().postevaluation + ' ' +
                          self.get_dir() + ' ' +
                          evaldir_args +
                          ' --challengedir ' +
                          os.path.join(chall.get_dir(),
                                       chall_dname) +
                          ' --stageprefix ' + self._eval_task_prefix_str +
                          ' --evaluationsuffix ' + self._eval_task_suffix_str)
            peval_name = os.path.basename(self.get_args().postevaluation)

            self.run_external_command(peval_name, cmd_to_run,
                                      True)
        else:
            self.set_error('Did not find any Evaluation tasks to summarize')

        # attempt to send postevaluation email
        try:
            self._emailer.send_postevaluation_email(self)
            self.append_to_email_log(self._emailer.get_message_log())
        except Exception as e:
            logger.exception('Caught exception trying to send '
                             'postevaluation email')
            self.append_to_email_log('Caught exception trying to send '
                                     'evaluation email ' + str(e) + '\n')

        # attempt to post evaluation results to website REST service
        try:
            self.post_target_stats_to_websiteservice(self.
                                                     generate_target_object())
        except Exception as e:
            logger.exception('Not a show stopper, but caught exception '
                             'trying to post results to website rest '
                             'service')
            self.append_to_email_log('\nCaught exception trying to post '
                                     'evaluation results to website rest '
                                     'service ' + str(e) + '\n')
        # assess the result
        self.end()
