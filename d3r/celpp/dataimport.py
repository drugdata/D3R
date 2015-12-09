# -*- coding: utf-8 -*-

__author__ = 'churas'

import logging
import os
from d3r.celpp.task import D3RTask
from d3r.celpp import util

logger = logging.getLogger(__name__)


class DataImportTask(D3RTask):
    """Represents DataImport Task
    This task downloads 3 files named new_release_structure_nonpolymer.tsv,
    new_release_structure_sequence.tsv, and
    new_release_crystallization_pH.tsv from the web
    """
    NONPOLYMER_TSV = "new_release_structure_nonpolymer.tsv"
    SEQUENCE_TSV = "new_release_structure_sequence.tsv"
    CRYSTALPH_TSV = "new_release_crystallization_pH.tsv"
    COMPINCHI_ICH = "Components-inchi.ich"

    def __init__(self, path, args):
        super(DataImportTask, self).__init__(path, args)
        self.set_name('dataimport')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._maxretries = 3
        self._retrysleep = 1

    def get_nonpolymer_tsv(self):
        """Returns path to new_release_structure_nonpolymer.tsv file
        :return: full path to DataImportTask.NONPOLYMER_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.NONPOLYMER_TSV)

    def get_sequence_tsv(self):
        """Returns path to new_release_structure_sequence.tsv file
        :return: full path to DataImportTask.SEQUENCE_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.SEQUENCE_TSV)

    def get_crystalph_tsv(self):
        """Returns path to new_release_crystallization_pH.tsv file
        :return: full path to DataImportTask.CRYSTALPH_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.CRYSTALPH_TSV)

    def get_components_inchi_file(self):
        return os.path.join(self.get_dir(),
                            DataImportTask.COMPINCHI_ICH)

    def can_run(self):
        """Determines if task can actually run

           The method verifies a `DataImportTask` does
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

           Downloads 3 files from url specified in `self._args.pdbfileurl`
           namely the ones mentioned in the constructor to get_dir() directory.
           Each download will be retried up to self._maxretries time
           which is set in constructor.  After which set_error()
           will be set with message if file is still unable to be downloaded.
           """
        super(DataImportTask, self).run()

        try:
            logger.debug('pdbfileurl set to ' + self._args.pdbfileurl)
        except AttributeError:
            self.set_error('cannot download files cause pdbfileurl not set')
            self.end()
            return

        try:
            logger.debug('compinchi set to ' + self._args.compinchi)
        except AttributeError:
            self.set_error('cannot download files cause compinchi not set')
            self.end()
            return

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        download_path = self.get_nonpolymer_tsv()
        url = self._args.pdbfileurl
        try:
            util.download_url_to_file(url +
                                      '/' + DataImportTask.NONPOLYMER_TSV,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            download_path = self.get_sequence_tsv()
            util.download_url_to_file(url +
                                      '/' + DataImportTask.SEQUENCE_TSV,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            download_path = self.get_crystalph_tsv()
            util.download_url_to_file(url +
                                      '/' + DataImportTask.CRYSTALPH_TSV,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            url = self._args.compinchi
            download_path = self.get_components_inchi_file()
            util.download_url_to_file(url +
                                      '/' + DataImportTask.COMPINCHI_ICH,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)
            self.end()
            return
        except Exception:
            logger.exception('Caught Exception trying to download file(s)')

        self.set_error('Unable to download file from ' +
                       url + ' to ' + download_path)

        # assess the result
        self.end()
