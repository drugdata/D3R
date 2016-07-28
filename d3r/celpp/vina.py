# -*- coding: utf-8 -*-

__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.chimeraprep import ChimeraProteinLigPrepTask

logger = logging.getLogger(__name__)


class AutoDockVinaTask(D3RTask):
    """Performs Auto Dock Vina docking

    """
    def __init__(self, path, args):
        super(AutoDockVinaTask, self).__init__(path, args)
        self.set_name('autodockvina')

        prep = ChimeraProteinLigPrepTask(path, args)
        self.set_stage(prep.get_stage() + 1)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain stderr/stdout files
           returns: list of files that can be uploaded.
        """
        # get the stderr/stdout files
        file_list = super(AutoDockVinaTask, self).get_uploadable_files()
        return file_list

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the `ProteinLigPrep` task
           has `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies a `AutoDockVinaTask` does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        # check blast
        ligprep = ChimeraProteinLigPrepTask(self._path, self._args)
        ligprep.update_status_from_filesystem()
        if ligprep.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + ligprep.get_name() + 'task' +
                        'has a status of ' + ligprep.get_status())
            self.set_error(ligprep.get_name() + ' task has ' +
                           ligprep.get_status() + ' status')
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
        """Runs AutoDockVinaTask after verifying chimeraproteinligprep was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes blastnfilter script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(AutoDockVinaTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('vina set to ' +
                         self.get_args().vina)
        except AttributeError:
            self.set_error('vina not set')
            self.end()
            return

        ligprep = ChimeraProteinLigPrepTask(self._path, self._args)

        #
        # vinadocking.py --structuredir <path to stage.3.chimeraprep> \
        # --outdir <path to stage.4.vina>
        #
        cmd_to_run = (self.get_args().vina + ' --structuredir ' +
                      ligprep.get_dir() +
                      ' --outdir ' + self.get_dir())

        vina_name = os.path.basename(self.get_args().vina)

        self.run_external_command(vina_name, cmd_to_run,
                                  True)
        # assess the result
        self.end()
