#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time

from d3r.celpp.task import D3RTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp import util
logger = logging.getLogger(__name__)

class ExternalDataSubmissionFactory(object):
    """Factory to create ExternalDataSubmissionObjects
    """


    def get_external_data_submissions(self):
        """Generate ExternalDataSubmission objects
           Method should examine the submission directory on the
           remote server via FtpFileTransfer(args.ftpconfig)
           for each directory under submission directory
           look for a celpp_weekXX_YYYY_dockedresults_ZZZZ.tar.gz file
           where ZZZZ matches name of directory under submission directory
           Regardless if its found or not create new ExternalDataSubmission object
           and pass ZZZZ for name and path to .tar.gz file as remotefile if
           the tar.gz file was found otherwise pass None
           along with args.  Append this object to list and
           return it
        pseudo code:
           filetransfer = self.get_file_transfer()
           dlist = filetransfer.list_dirs(ft.get_ftp_remote_submission_dir()
           for d in dlist:
             for x in filetransfer.list_files:
               if x starts with(chall_dir_name) and ends in .tar.gz
                 # create new external submission object and add to list
           return #list
        """
        return []

class ExternalDataSubmissionTask(D3RTask):
    """Downloads external user docking Submissions
    """

    def __init__(self, path, name, remotefile, args):
        super(ExternalDataSubmissionTask, self).__init__(path, args)
        self.set_name(name)
        self.set_stage(EvaluationTaskFactory.DOCKSTAGE)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self.set_remote_challenge_data_package(remotefile)

    def set_remote_challenge_data_package(self, remotefile):
        self._remote_challenge_file = remotefile

    def can_run(self):
        """Determines if task can actually run
           by seeing if this task has already been run
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

    def _untar_challenge_data_package(self):
        pass

    def _get_summary_of_docked_results(self):
        pass

    def _download_external_submission(self):
        """Downloads external submissions
           pseudo code:
           ft.download_file(self.get_remote_challenge_data_package())
        """
        pass

    def run(self):
        """Runs ExternalDataSubmission task

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(ExternalDataSubmissionTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return
        try:
            # self._download_remote_challenge_data_package()
            # self._untar_challenge_data_package()
            # summary = self._get_summary_of_docked_results()
            # self.append_to_email_log(summary)
        except Exception as e:
            logger.exception('Caught exception')
            self.set_error('Caught exception ' + str(e))

        # assess the result
        self.end()
