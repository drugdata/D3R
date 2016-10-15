__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.challengedata import ChallengeDataTask

logger = logging.getLogger(__name__)


class ChimeraProteinLigPrepTask(D3RTask):
    """Performs preparation of PDB and inchi files

    """

    FINAL_LOG = 'final.log'

    def __init__(self, path, args):
        super(ChimeraProteinLigPrepTask, self).__init__(path, args)
        self.set_name('chimeraprep')

        # Make stage number one higher then BlastNFilter Stage
        chall = ChallengeDataTask(path, args)
        self.set_stage(chall.get_stage() + 1)

        self.set_status(D3RTask.UNKNOWN_STATUS)

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files in addition to stderr/stdout files
           if they are found on the filesystem

           final.log

           :returns: list of files that can be uploaded.
        """
        # get the stderr/stdout files
        file_list = super(ChimeraProteinLigPrepTask,
                          self).get_uploadable_files()

        out_dir = self.get_dir()

        final_log = os.path.join(out_dir, ChimeraProteinLigPrepTask.FINAL_LOG)
        if os.path.isfile(final_log):
            file_list.append(final_log)

        return file_list

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the `ChallengeDataTask` task
           has `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies a `ProteinLigPrepTask` does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        # check blast
        chall = ChallengeDataTask(self._path, self._args)
        chall.update_status_from_filesystem()
        if chall.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + chall.get_name() + 'task' +
                        'has a status of ' + chall.get_status())
            self.set_error(chall.get_name() + ' task has ' +
                           chall.get_status() + ' status')
            return False

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
        """Runs ChimeraProteinLigPrepTask after verifying BlastNFilterTask
           was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes blastnfilter script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(ChimeraProteinLigPrepTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('chimeraprep set to ' +
                         self.get_args().chimeraprep)
        except AttributeError:
            self.set_error('chimeraprep not set')
            self.end()
            return

        try:
            logger.debug('pdbdb set to ' +
                         self.get_args().pdbdb)
        except AttributeError:
            self.set_error('pdbdb not set')
            self.end()
            return

        rdkitpython = ''
        try:
            logger.debug('rdkitpython set to  ' +
                         self.get_args().rdkitpython)
            rdkitpython = self.get_args().rdkitpython
        except AttributeError:
            logger.debug

        chall = ChallengeDataTask(self._path, self._args)

        #
        # chimeraprep.py --candidatedir <path to stage.3.challengedata> \
        # --outdir <path to stage.4.chimeraligprep>

        challdir = os.path.join(chall.get_dir(),
                                chall.get_celpp_challenge_data_dir_name())

        cmd_to_run = (self.get_args().chimeraprep + ' --candidatedir ' +
                      challdir +
                      ' --rdkitpython \'' + rdkitpython + '\'' +
                      ' --pdbdb ' + self.get_args().pdbdb +
                      ' --outdir ' + self.get_dir())

        chimeraprep_name = os.path.basename(self.get_args().chimeraprep)

        self.run_external_command(chimeraprep_name, cmd_to_run,
                                  True)
        # assess the result
        self.end()
