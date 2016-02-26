__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask

logger = logging.getLogger(__name__)


class ProteinLigPrepTask(D3RTask):
    """Performs preparation of PDB and inchi files

    """

    def __init__(self, path, args):
        super(ProteinLigPrepTask, self).__init__(path, args)
        self.set_name('proteinligprep')
        self.set_stage(3)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the `BlastNFilterTask` task
           has `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies a `ProteinLigPrepTask` does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        # check blast
        blastnfilter = BlastNFilterTask(self._path, self._args)
        blastnfilter.update_status_from_filesystem()
        if blastnfilter.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + blastnfilter.get_name() + 'task' +
                        'has a status of ' + blastnfilter.get_status())
            self.set_error(blastnfilter.get_name() + ' task has ' +
                           blastnfilter.get_status() + ' status')
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
        """Runs ProteinLigPrepTask after verifying blastnfilter was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes blastnfilter script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(ProteinLigPrepTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('proteinligprep set to ' +
                         self.get_args().proteinligprep)
        except AttributeError:
            self.set_error('proteinligprep not set')
            self.end()
            return

        try:
            logger.debug('pdbdb set to ' +
                         self.get_args().pdbdb)
        except AttributeError:
            self.set_error('pdbdb not set')
            self.end()
            return

        blastnfilter = BlastNFilterTask(self._path, self._args)

        #
        # proteinligprep.py --candidatedir <path to stage.2.blastnfilter> \
        # --outdir <path to stage.3.proteinligprep>
        cmd_to_run = (self.get_args().proteinligprep + ' --candidatedir ' +
                      blastnfilter.get_dir() +
                      ' --pdbdb ' + self.get_args().pdbdb +
                      ' --outdir ' + self.get_dir())

        proteinligprep_name = os.path.basename(self.get_args().proteinligprep)

        self.run_external_command(proteinligprep_name, cmd_to_run,
                                  True)
        # assess the result
        self.end()
