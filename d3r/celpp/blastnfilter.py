# -*- coding: utf-8 -*-

__author__ = 'churas'

import os
import logging
import re

from d3r.celpp import util
from d3r.celpp.task import D3RTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from d3r.celpp.dataimport import DataImportTask


logger = logging.getLogger(__name__)


class BlastNFilterSummary:
    """Represents summary of a BlastNFilterTask invocation
    """
    def __init__(self, path):

        if path is None:
            logger.error('Path passed into constructor is None type')
            self._path = ''
        else:
            self._path = path

        self._complexes = 0
        self._dockable_complexes = 0
        self._dockable_monomers = 0
        self._targets_found = 0
        self._parse_summary_file()

    def _parse_summary_file(self):
        """Parses summary file to obtain needed fields in object

           WARNING: This method assumes a well formed summary.txt file and
           if its not then output maybe wonky
        """
        if not os.path.isdir(self._path):
            logger.error(self._path + ' is not a directory')
            return

        # open summary.txt file and extract
        # complexes, dockable complexes, dockable monomers and Targets found
        # below is example summary.txt file with parts omitted
        #
        # INPUT SUMMARY
        #   entries:                             135
        #   complexes:                            83
        #   dockable complexes:                   46
        #   monomers:                             90
        #   dockable monomers:                    38
        #   multimers:                            45
        #   dockable multimers:                    8
        # .
        # . <THERE IS OTHER DATA HERE>
        # .
        # OUTPUT SUMMARY
        #   Targets found:                        33
        # .
        # .
        summary_txt = os.path.join(self._path, BlastNFilterTask.SUMMARY_TXT)
        if not os.path.isfile(summary_txt):
            logger.error('No ' + summary_txt + ' file found')
            return

        complexes_pat = re.compile("^ +complexes: +")
        dock_complex_pat = re.compile("^ +dockable complexes: +")
        dock_mon_pat = re.compile("^ +dockable monomers: +")
        targets_pat = re.compile("^ +Targets found: +")

        logger.debug('Parsing file ' + summary_txt)

        f = open(summary_txt, 'r')
        for line in f:
            if re.match(complexes_pat, line):
                self.set_complexes(re.sub("^.*: +", "", line.rstrip()))
            elif re.match(dock_complex_pat, line):
                self.set_dockable_complexes(re.sub("^.*: +", "",
                                                   line.rstrip()))
            elif re.match(dock_mon_pat, line):
                self.set_dockable_monomers(re.sub("^.*: +", "", line.rstrip()))
            elif re.match(targets_pat, line):
                self.set_targets_found(re.sub("^.*: +", "", line.rstrip()))
        f.close()

    def get_complexes(self):
        return self._complexes

    def set_complexes(self, val):
        self._complexes = val

    def get_dockable_complexes(self):
        return self._dockable_complexes

    def set_dockable_complexes(self, val):
        self._dockable_complexes = val

    def get_dockable_monomers(self):
        return self._dockable_monomers

    def set_dockable_monomers(self, val):
        self._dockable_monomers = val

    def get_targets_found(self):
        return self._targets_found

    def set_targets_found(self, val):
        self._targets_found = val

    def get_week_number(self):
        """Parses path for week number
           :returns: week number as string
        """
        return util.get_celpp_week_number_from_path(
            os.path.dirname(self._path))

    def get_year(self):
        """Parses year from path passed into constructor
           :returns: year as string
        """
        return util.get_celpp_year_from_path(self._path)

    def get_csv(self):
        """Returns comma separated string of values parsed from summary.txt

           :returns: Week #, Year, # complexes, # dockable complexes,
                   # dockable monomers, # targets found
        """
        return (self.get_week_number() + ',' + self.get_year() + ',' +
                str(self.get_complexes()) + ',' +
                str(self.get_dockable_complexes()) +
                ',' + str(self.get_dockable_monomers()) + ',' +
                str(self.get_targets_found()))


class BlastNFilterTask(D3RTask):
    """Performs blast and filter of sequences

    """
    SUMMARY_TXT = "summary.txt"
    DOCKABLE_XSLX = "dockable.xlsx"
    BLASTNFILTER_LOG = "blastnfilter.log"
    TASK_NAME = 'blastnfilter'

    def __init__(self, path, args):
        super(BlastNFilterTask, self).__init__(path, args)
        self.set_name(BlastNFilterTask.TASK_NAME)

        # set the stage to be 1 higher then Data Import Stage
        dataimport = DataImportTask(path, args)
        self.set_stage(dataimport.get_stage() + 1)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def get_blastnfilter_summary(self):
        """Factory method that Generates BlastNFilterSummary object

           by parsing data from this object
        """
        return BlastNFilterSummary(self.get_dir())

    def get_blastnfilter_summary_file(self):
        """Returns path to summary.txt file
        """
        return os.path.join(self.get_dir(), BlastNFilterTask.SUMMARY_TXT)

    def _parse_blastnfilter_output_for_hit_stats(self):
        """Examines output directory of blastnfilter.py for stats on run

           This method looks at the output directory and counts the number
           of .csv files found.  Each .csv file corresponds to a target.
           Within each .csv file the entries under Test Information correspond
           to candidates.  This method will output number of candidates for
           that target by calling get_candidate_count_from_csv method
        """

        txt_list = self.get_txt_files()
        summary_out = '\nOutput from ' + BlastNFilterTask.SUMMARY_TXT + '\n'
        summary_txt = os.path.join(self.get_dir(),
                                   BlastNFilterTask.SUMMARY_TXT)
        try:
            if os.path.isfile(summary_txt):
                f = open(summary_txt, 'r')
                for entry in f:
                    summary_out = (summary_out + entry)
                f.close()
        except:
            logger.warning('Error reading ' + summary_txt + ' file')

        return '\n# txt files found: ' + str(len(txt_list)) +\
               '\n' + summary_out

    def get_txt_files(self, addfullpath=False):
        """ Gets txt files in task directory (just the names) skiping summary.txt
        :return:list of txt file names
        """
        out_dir = self.get_dir()
        txt_list = []
        try:
            for entry in os.listdir(out_dir):
                if entry == BlastNFilterTask.SUMMARY_TXT:
                    continue

                if entry.endswith('.txt'):
                    if addfullpath is True:
                        txt_list.append(os.path.join(out_dir, entry))
                    else:
                        txt_list.append(entry)
        except OSError:
            logger.warning('Caught exception trying to look for .txt files')

        return txt_list

    def get_uploadable_files(self):
        """Returns a list of files that can be uploaded to remote server

           Provides a way for a D3RTask to define what files should
           uploaded to remote server by any D3RTaskUploader.  The default
           implementation returns an empty list
           :returns: list of files that can be uploaded.  The files should
           have full paths
        """

        # get the stderr/stdout files
        file_list = super(BlastNFilterTask, self).get_uploadable_files()

        out_dir = self.get_dir()
        # add txt files
        file_list.extend(self.get_txt_files(True))

        summary_file = os.path.join(out_dir,
                                    BlastNFilterTask.SUMMARY_TXT)
        if os.path.isfile(summary_file):
            file_list.append(summary_file)

        dockable_file = os.path.join(out_dir,
                                     BlastNFilterTask.DOCKABLE_XSLX)

        if os.path.isfile(dockable_file):
            file_list.append(dockable_file)

        bnf_log_file = os.path.join(out_dir, BlastNFilterTask.BLASTNFILTER_LOG)

        if os.path.isfile(bnf_log_file):
            file_list.append(bnf_log_file)
        return file_list

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

        try:
            loglevel = self.get_args().loglevel
        except AttributeError:
            logger.debug('No log level set in arguments using WARNING')
            loglevel = 'WARNING'

        # verify sequence.tsv file exists on filesystem.
        # if not fall back to oldsequence.tsv file
        sequencetsv = data_import.get_sequence_tsv()
        if not os.path.isfile(sequencetsv):
            logger.warning(sequencetsv + ' file not found. falling '
                                         'back to old file')
            self.append_to_email_log('\n ' + sequencetsv + ' file not found ' +
                                     'falling back to ' +
                                     data_import.get_oldsequence_tsv() + '\n')
            sequencetsv = data_import.get_oldsequence_tsv()

        cmd_to_run = (self.get_args().blastnfilter + ' --nonpolymertsv ' +
                      data_import.get_nonpolymer_tsv() +
                      ' --sequencetsv ' +
                      sequencetsv +
                      ' --pdbblastdb ' +
                      make_blastdb.get_dir() +
                      ' --compinchi ' +
                      data_import.get_components_inchi_file() +
                      ' --crystalpH ' +
                      data_import.get_crystalph_tsv() +
                      ' --pdbdb ' +
                      self.get_args().pdbdb +
                      ' --log ' +
                      loglevel +
                      ' --outdir ' + self.get_dir())

        blastnfilter_name = os.path.basename(self.get_args().blastnfilter)

        self.run_external_command(blastnfilter_name,
                                  cmd_to_run, False,)

        self.set_status(D3RTask.COMPLETE_STATUS)

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
