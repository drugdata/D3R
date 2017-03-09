__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.task import D3RParameters
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp import util


logger = logging.getLogger(__name__)


class PathNotDirectoryError(Exception):
    """Path is not a directory Error
    """
    pass


class PostEvaluationTask(D3RTask):
    """Generates summary of all evaluations run
    """
    RMSD_TXT = 'RMSD.txt'
    RMSD_PICKLE = 'RMSD.pickle'
    TASK_NAME = 'postevaluation'
    FINAL_LOG = 'final.log'
    CSV_SUFFIX = '.csv'

    def __init__(self, path, args):
        super(PostEvaluationTask, self).__init__(path, args)
        self.set_name(PostEvaluationTask.TASK_NAME)
        et = EvaluationTask('/foo', '', None, None)
        self._eval_task_prefix_str = et.get_dir_name()
        self.set_stage(et.get_stage() + 1)
        self._emailer = None
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def set_evaluation_emailer(self, emailer):
        self._emailer = emailer

    def get_all_completed_evaluation_tasks(self):
        """Look for all completed evaluation tasks
        """
        task_list = []
        prefix_len = len(self._eval_task_prefix_str)
        base_dir = self.get_path()
        for entry in os.listdir(base_dir):
            if not entry.endswith(EvaluationTaskFactory.SCORING_SUFFIX):
                continue
            full_path = os.path.join(base_dir, entry)
            if not os.path.isdir(full_path):
                logger.info('Suffix matched, but is not directory: ' +
                            full_path)
                continue
            task = EvaluationTask(base_dir, entry[prefix_len:], None,
                                  self.get_args())
            task.update_status_from_filesystem()
            if task.get_status() != D3RTask.COMPLETE_STATUS:
                logger.info('Skipping ' + task.get_name() + ' task ' +
                            'because ' + task.get_name() + ' task' +
                            ' has a status of ' + task.get_status())
                continue
            task_list.append(task)
        return task_list

    def _get_evaluationdir_args(self):
        """Generates a string containing all the evaluation
           tasks to pass to post_evaluation.py script
        :returns: string of format --evaluationdir <path1> --evaluationdir <path2>
        """
        return ''

    def get_all_csv_files_in_dir(self):
        """Gets all the files ending with CSV_SUFFIX
           in task directory
        """
        csv_list = []
        out_dir = self.get_dir()
        for entry in os.listdir(out_dir):
            if entry.endswith(PostEvaluationTask.CSV_SUFFIX):
                csv_list.append(os.path.join(out_dir, entry))
        return csv_list

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
        # --challengedir <path to challenge dir>
        # --outdir <path to stage.5.glide.evaluation>
        #
        cmd_to_run = (self.get_args().postevaluation + ' ' +
                      self._get_evaluationdir_args() +
                      ' --outdir ' + self.get_dir())

        peval_name = os.path.basename(self.get_args().postevaluation)

        self.run_external_command(peval_name, cmd_to_run,
                                  True)

        # attempt to send postevaluation email
        # try:
        #    # self._emailer.send_evaluation_email(self)
        #    # self.append_to_email_log(self._emailer.get_message_log())
        # except Exception as e:
        #    logger.exception('Caught exception trying to send evaluation '
        #                     'email')
        #    self.append_to_email_log('Caught exception trying to send '
        #                             'evaluation email ' + str(e) + '\n')

        # assess the result
        self.end()
