#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import tarfile
import shutil
import time

from d3r.celpp.task import D3RTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.filetransfer import FtpFileTransfer

logger = logging.getLogger(__name__)


class ChallengePackageDownloadError(Exception):
    """Raised when there was an error downloading
       Challenge data package
    """
    pass


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
        :raises AttributeError: if `get_file_transfer()` is None
        """
        dlist = self._file_transfer.list_dirs(remote_dir)
        if dlist is None:
            logger.debug('No directories returned')
            return []
        logger.debug('Found ' + str(len(dlist)) + ' directories on ' +
                     remote_dir)
        return dlist

    def _get_challenge_data_package_file(self, remote_dir, dir_name):
        """Gets challenge data package file under `remote_dir` / `dir_name`
           if it exists
           :returns: Path to remote challenge file upon success or None if
           not found
           :raises AttributeError if `get_file_transfer()` is None
        """
        ft = self.get_file_transfer()
        flist = ft.list_files(os.path.normpath(os.path.join(remote_dir,
                                                            dir_name)))
        if flist is None:
            logger.info('No files found in ' + dir_name)
            return None
        logger.debug('Found ' + str(len(flist)) + ' files in ' + dir_name)
        chall_fname = self._get_challenge_package_results_file_name(dir_name)
        logger.debug('Looking for ' + chall_fname + ' file in directory')

        for entry in flist:
            if entry == chall_fname:
                logger.info('Found matching entry ' + entry)
                return os.path.normpath(os.path.join(remote_dir,
                                                     dir_name,
                                                     chall_fname))
            else:
                logger.debug('Encountered non challenge file: ' + entry)

        return None

    def _remove_latest_txt(self):
        """Removes the latest.txt file from ftp server if found
        """
        try:
            ft = self.get_file_transfer()
            latest_txt = os.path.join(ft.get_remote_challenge_dir(),
                                      ChallengeDataTask.LATEST_TXT)
            logger.info('Attempting to remove ' + latest_txt)
            val = ft.delete_file(latest_txt)
            logger.info('Return value from delete call ' + str(val))
        except Exception:
            logger.exception('Caught exception trying to remove latest.txt')

    def get_external_data_submissions(self):
        """Generate ExternalDataSubmission objects
           Method should examine the submission directory on the
           remote server via FtpFileTransfer(args.ftpconfig)
           for each directory under submission directory
           look for a celpp_weekXX_YYYY_dockedresults_ZZZZ.tar.gz file
           where ZZZZ matches name of directory under submission directory
           Regardless if its found or not create new ExternalDataSubmission
           object and pass ZZZZ for name and path to .tar.gz file as
           remotefile if the tar.gz file was found otherwise pass None
           along with args.  Append this object to list and
           return it
        """
        task_list = []
        try:
            self._file_transfer.connect()
            self._remove_latest_txt()
            subdir = self._file_transfer.get_remote_submission_dir()
            dlist = self._get_submission_dirs(subdir)
            for d in dlist:
                chall_file = self._get_challenge_data_package_file(subdir, d)

                if chall_file is not None:
                    et = ExternalDataSubmissionTask(self.get_path(), d,
                                                    chall_file,
                                                    self.get_args())
                    logger.info('Added ExternalData Submission Task: ' +
                                et.get_dir_name())
                    task_list.append(et)
        except Exception:
            logger.exception('Caught exception')
        finally:
            try:
                self._file_transfer.disconnect()
            except Exception:
                logger.exception('Caught exception disconnecting')

        return task_list


class ExternalDataSubmissionTask(D3RTask):
    """Downloads external user docking Submissions
    """

    def __init__(self, path, name, remotefile, args):
        super(ExternalDataSubmissionTask, self).__init__(path, args)
        self.set_name(name + EvaluationTask.EXT_SUBMISSION_SUFFIX)
        self.set_stage(EvaluationTaskFactory.DOCKSTAGE)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self.set_remote_challenge_data_package(remotefile)
        self._maxretries = 3
        self._retrysleep = 30

    def set_download_max_retry_count(self, numretries):
        """Sets number of retries to perform before giving up on
           download
           :param numretries: int representing number of retries
        """
        self._maxretries = int(numretries)

    def set_download_retry_sleep(self, retrysleep):
        """Sets number of seconds to wait before retrying download
           :param retrysleep: int containing number of seconds
        """
        self._retrysleep = int(retrysleep)

    def set_remote_challenge_data_package(self, remotefile):
        self._remote_challenge_file = remotefile

    def get_remote_challenge_data_package(self):
        return self._remote_challenge_file

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

    def _is_tarmembername_safe(self, tname, chall_name):
        if tname.startswith(os.sep):
            msg = ('Skipping, path starts with / : ' +
                   tname)
            logger.warning(msg)
            self.append_to_email_log(msg + '\n')
            return False

        if not tname.startswith(chall_name):
            msg = ('Skipping, path does not conform: ' + tname)
            logger.warning(msg)
            self.append_to_email_log(msg + '\n')
            return False

        if tname.find('..') > -1:
            msg = ('Skipping, found .. in path: ' + tname)
            logger.warning(msg)
            self.append_to_email_log(msg + '\n')
            return False
        return True

    def _untar_challenge_data_package(self, chall_file):
        """Untars tar file into task working directory
        """
        try:
            tar = tarfile.open(os.path.join(self.get_dir(), chall_file))
            chall_name = chall_file.replace(ChallengeDataTask.TAR_GZ_SUFFIX,
                                            '')
            for t in tar.getmembers():
                if self._is_tarmembername_safe(t.name, chall_name) is False:
                    continue

                prefix_dir_removed = t.name.replace(chall_name + '/', '', 1)
                logger.debug('Found entry in tar ' + t.name)
                if t.isdir():
                    logger.debug('Found directory ' + prefix_dir_removed)
                    continue
                if t.isfile():
                    logger.debug('Extracting file ' + prefix_dir_removed)
                    tar.extract(t, path=self.get_dir())
                    continue
                msg = ('Ignoring non dir/file entry in tar ' +
                       prefix_dir_removed)
                logger.warning(msg)
                self.append_to_email_log(msg + '\n')
            return chall_name
        finally:
            try:
                tar.close()
            except Exception:
                logger.exception('Caught exception trying to close tar')

    def _move_challenge_data_package_into_task_dir(self, chall_name):
        """Moves the contents of uncompressed challenge data package into
           `get_dir()` of task
        """
        package_dir = os.path.join(self.get_dir(), chall_name)
        for entry in os.listdir(package_dir):
            fullpath = os.path.join(package_dir, entry)
            logger.debug('Moving ' + fullpath)
            shutil.move(fullpath, os.path.join(self.get_dir(), entry))
        os.rmdir(package_dir)

    def _get_summary_of_docked_results(self):
        """Right now just returns an empty string
        """
        return ''

    def _download_remote_challenge_data_package_with_retry(self):
        """With retries attempts to download remote challenge data
           package
           :returns: string containing path to remote challenge data
                     package file or None upon failure.
        """
        chall_name = os.path.basename(self.get_remote_challenge_data_package())
        count = 0
        while count <= self._maxretries:
            logger.debug('Try # ' + str(count) + ' of ' +
                         str(self._maxretries) + ' to download ' +
                         chall_name)
            if count > 0:
                self.append_to_email_log('Try # ' + str(count) + ' of ' +
                                         str(self._maxretries) +
                                         ' to download ' + chall_name + '\n')

            try:
                res = self._download_remote_challenge_data_package(chall_name)
                if res is True:
                    logger.debug('Successfully downloaded')
                    return chall_name
                logger.debug('Download failed, sleeping ' +
                             str(self._retrysleep) + ' seconds')
                time.sleep(self._retrysleep)
            except Exception:
                logger.exception('Caught exception, retrying download')
            finally:
                count += 1

        raise ChallengePackageDownloadError('Unable to download ' + chall_name)

    def _download_remote_challenge_data_package(self, chall_name):
        """Downloads external submissions
           pseudo code:
           ft.download_file(self.get_remote_challenge_data_package())
        """
        try:

            ft = self.get_file_transfer()
            ft.connect()
            localfile = os.path.join(self.get_dir(), chall_name)
            logger.debug('Downloading ' +
                         self.get_remote_challenge_data_package() + ' to ' +
                         localfile)
            return ft.download_file(self.get_remote_challenge_data_package(),
                                    localfile)
        finally:
            try:
                ft.disconnect()
            except Exception:
                logger.exception('Caught exception trying to disconnect')

    def _runtask(self):
        """Performs the download and untar of challenge data package
        """
        chall_file = self._download_remote_challenge_data_package_with_retry()
        chall_name = self._untar_challenge_data_package(chall_file)
        self._move_challenge_data_package_into_task_dir(chall_name)
        summary = self._get_summary_of_docked_results()
        self.append_to_email_log(summary)

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
            self._runtask()
        except Exception as e:
            logger.exception('Caught exception')
            self.set_error('Caught exception ' + str(e))

        # assess the result
        self.end()
