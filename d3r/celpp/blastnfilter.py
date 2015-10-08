# -*- coding: utf-8 -*-

__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp.dataimport import DataImportTask

logger = logging.getLogger(__name__)


class BlastNFilterTask(D3RTask):
    """Performs blast and filter of sequences

    """

    def __init__(self, path, args):
        super(BlastNFilterTask, self).__init__(path, args)
        self.set_name('blastnfilter')
        self.set_stage(2)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def _parse_blastnfilter_output_for_hit_stats(self):
        """Examines output directory of blastnfilter.py for stats on run

           This method looks at the output directory and counts the number
           of .csv files found.  Each .csv file corresponds to a target.
           Within each .csv file the entries under Test Information correspond
           to candidates.  This method will output number of candidates for
           that target by calling get_candidate_count_from_csv method
        """

        txt_list = self.get_txt_files()
        target_list = '\n\nTarget\n'
        for entry in txt_list:
            target_list = (target_list + entry.replace('.txt', '') + '\n')

        return '\n# targets found: ' + str(len(txt_list)) + '\n' + target_list

    def get_txt_files(self):
        """ Gets txt files in task directory (just the names) skiping summary.txt
        :return:list of txt file names
        """
        out_dir = self.get_dir()
        txt_list = []
        for entry in os.listdir(out_dir):
            if entry == 'summary.txt':
                continue

            if entry.endswith('.txt'):
                txt_list.append(entry)

        return txt_list

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the `MakeBlastDBTask` and
           `DataImportTask` task have
           `D3RTask.COMPLETE_STATUS` for status.  The method then
           verifies a `BlastNFilterTask` does not already exist.
             If above is not true then self.set_error() is set
             with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        # check blast
        make_blastdb = MakeBlastDBTask(self._path, self._args)
        make_blastdb.update_status_from_filesystem()
        if make_blastdb.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + ' task ' +
                        'because ' + make_blastdb.get_name() + ' task' +
                        'has a status of ' + make_blastdb.get_status())
            self.set_error(make_blastdb.get_name() + ' task has ' +
                           make_blastdb.get_status() + ' status')
            return False

        # check data import
        data_import = DataImportTask(self._path, self._args)
        data_import.update_status_from_filesystem()
        if data_import.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + ' task ' +
                        'because ' + data_import.get_name() + ' task' +
                        'has a status of ' + data_import.get_status())
            self.set_error(data_import.get_name() + ' task has ' +
                           data_import.get_status() + ' status')
            return False

        # check blast is not complete and does not exist

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
        """Runs blastnfilter task after verifying dataimport was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes blastnfilter script and
           postanalysis script.  Upon completion results are
           analyzed and success or error status is set
           appropriately and D3RTask.end is invoked
           """
        super(BlastNFilterTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        data_import = DataImportTask(self._path, self._args)

        make_blastdb = MakeBlastDBTask(self._path, self._args)

        cmd_to_run = (self.get_args().blastnfilter + ' --nonpolymertsv ' +
                      data_import.get_nonpolymer_tsv() +
                      ' --sequencetsv ' +
                      data_import.get_sequence_tsv() +
                      ' --pdbblastdb ' +
                      make_blastdb.get_dir() +
                      ' --compinchi ' +
                      data_import.get_components_inchi_file() +
                      ' --outdir ' + self.get_dir())

        blastnfilter_name = os.path.basename(self.get_args().blastnfilter)

        self.run_external_command(blastnfilter_name,
                                  cmd_to_run, True)

        cmd_to_run = (self.get_args().postanalysis + ' --compinchi ' +
                      data_import.get_components_inchi_file() + ' ' +
                      self.get_dir())

        postanalysis_name = os.path.basename(self.get_args().postanalysis)

        self.run_external_command(postanalysis_name,
                                  cmd_to_run, False)

        try:
            # examine output to get candidate hit count DR-12
            hit_stats = self._parse_blastnfilter_output_for_hit_stats()
            if hit_stats is not None:
                self.append_to_email_log(hit_stats)
        except Exception:
            logger.exception("Error caught exception")

        # assess the result
        self.end()
