#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time

from d3r.celpp.task import D3RTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.filetransfer import FtpFileTransfer
from d3r.celpp import util

logger = logging.getLogger(__name__)

class ExternalDataSubmissionFactory(object):
    """Factory to create ExternalDataSubmissionObjects
    """

    DOCKEDRESULTS = '_dockedresults_'

    def __init__(self, path, args):
        """Constructor
        """
        try:
            logger.debug('ftpconfig set to ' + args.ftpconfig)
            self._file_transfer = FtpFileTransfer(args.ftpconfig)
        except Exception:
            logger.exception('Caught exception')
            self._file_transfer = None
        self._path = path
        self._args = args
        ctask = ChallengeDataTask(path, args)
        self._chall_dir_name = ctask.get_celpp_challenge_data_dir_name()

    def get_args(self):
        return self._args

    def get_path(self):
        return self._path

    def _get_challenge_package_results_file_name(self, dir_name):
        return (self._chall_dir_name +
                ExternalDataSubmissionFactory.DOCKEDRESULTS +
                dir_name + ChallengeDataTask.TAR_GZ_SUFFIX)

    def set_file_transfer(self, filetransfer):
        """Sets file transfer
        """
        self._file_transfer = filetransfer

    def get_file_transfer(self):
        """Gets file transfer
        """
        return self._file_transfer

    def _get_submission_dirs(self, remote_dir):
        """Gets list of directories under `remote_dir`
        :param remote_dir: path on remote server to examine
        :returns: list of directory names without path prefix
        """
        try:
            dlist = self._file_transfer.list_dirs(remote_dir)
            if dlist is None:
                logger.debug('No directories returned')
                return []
            logger.debug('Found ' + str(len(dlist)) + ' directories on ' +
                         remote_dir)
            return dlist
        except Exception:
            logger.exception('Caught exception')
        return []

    def _get_challenge_data_package_file(self, remote_dir, dir_name):
        ft = self.get_file_transfer()
        flist = ft.list_files(os.path.normpath(remote_dir, dir_name))
        if flist is None:
            logger.info('No files found in ' + dir_name)
            return None
        logger.debug('Found ' + str(len(flist)) + ' files in ' + dir_name)
        chall_fname = self._get_challenge_package_results_file_name(dir_name)
        logger.debug('Looking for ' + chall_fname + ' file in directory')

        for entry in flist:
            if entry == chall_fname:
                logger.debug('Found matching entry')
                return os.path.normpath(remote_dir, chall_fname)

        # if no matching files were found should we mention any ones
        # that have a similar name or that have been uploaded since
        # start of celpp week?

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
        """
        subdir = self._file_transfer.get_ftp_remote_submission_dir()
        dlist = self._get_submission_dirs(subdir)
        task_list = []
        for d in dlist:
            chall_file = self._get_challenge_data_package_file(subdir, d)

            if chall_file is not None:
                et = ExternalDataSubmissionTask(self.get_path(), d,
                                                chall_file, self.get_args())
                logger.info('Added ExternalData Submission Task: ' + et.get_dir_name())
                task_list.append(et)

        return task_list

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
        return ''

    def _download_remote_challenge_data_package(self):
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

            self._download_remote_challenge_data_package()
            self._untar_challenge_data_package()
            summary = self._get_summary_of_docked_results()
            self.append_to_email_log(summary)
        except Exception as e:
            logger.exception('Caught exception')
            self.set_error('Caught exception ' + str(e))

        # assess the result
        self.end()
