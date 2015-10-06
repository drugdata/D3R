# -*- coding: utf-8 -*-

__author__ = 'churas'

import logging
import urllib
import os
from d3r.celpp.task import D3RTask

logger = logging.getLogger(__name__)


class CompInchiDownloadTask(D3RTask):
    """Downloads Components-inchi.inchi file
    """

    def __init__(self, path, args):
        super(CompInchiDownloadTask, self).__init__(path, args)
        self.set_name('compinchi')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._maxretries = 3

    def get_components_inchi_file(self):
        return os.path.join(self.get_dir(), 'Components-inchi.ich')

    def can_run(self):
        """Determines if task can actually run

           The method verifies a `CompInchiDownloadTask` does
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
        """Downloads Components-inchi.ich

           Downloads data from url specified in self._args.compinchi to
           get_components_inchi_file()/Components-inchi.ich file.  Download
           will be retried up to self._maxretries time which is set in
           constructor.  After which set_error() will be set with message
           if file is still unable to be downloaded.
           """
        super(CompInchiDownloadTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return
        download_path = self.get_components_inchi_file()
        count = 1
        while count <= self._maxretries:
            logger.debug('Try # ' + str(count) + ' of ' +
                         str(self._maxretries) + ' to download ' +
                         download_path + ' from ' + self._args.compinchi)
            try:
                (f, info) = urllib.urlretrieve(self._args.compinchi,
                                               filename=download_path)
                self.end()
                return
            except Exception:
                logger.exception('Caught Exception trying to download file')
            count += 1

        self.set_error('Unable to download file from ' +
                       self._args.compinchi)

        # assess the result
        self.end()
