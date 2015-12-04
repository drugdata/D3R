__author__ = 'churas'

import os
import logging

from d3r.celpp.task import D3RTask
from d3r.celpp.proteinligprep import ProteinLigPrepTask

logger = logging.getLogger(__name__)

class PathNotDirectoryError(Exception):
    """Path is not a directory Error
    """
    pass

class ScoringTaskFactory(object):
    """Factory class to generate ScoringTask objects

       This factory examines a celpp week directory for
       all docking tasks.  The code then generates
       ScoringTask objects for all eligible docking tasks
    """
    DOCKSTAGE = 4
    STAGE_FOUR_PREFIX = D3RTask.STAGE_DIRNAME_PREFIX + '.' +\
                        str(DOCKSTAGE) + '.'
    SCORING_SUFFIX = 'scoring'
    WEB_DATA_SUFFIX = 'webdata'

    def __init__(self, path, theargs):
        """Constructor
        """
        self.set_path(path)
        self.set_args(theargs)

    def set_args(self, theargs):
        """ Sets args
        :param theargs: arguments to set
        """
        self._args = theargs

    def get_args(self):
        """Gets args passed into constructor or via set_args()
        """
        return self._args

    def set_path(self, path):
        """Sets path used to look for docking tasks
        """
        self._path = path

    def get_path(self):
        """Gets path used to look for docking tasks
        """
        return self._path

    def get_scoring_tasks(self):
        """Generate ScoringTasks

           This method examines the path directory
           set via the constructor or set_path() method
           for all stage 4 tasks excluding tasks
           that end with 'webdata'  A ScoringTask
           object is created for each of these tasks
           and returned in a list.
           :return: list of ScoringTask objects or empty list if none found
        """
        path = self.get_path()
        logger.debug('Examining ' + path + ' for docking tasks')
        if not os.path.isdir(path):
            raise PathNotDirectoryError(path + ' is not a directory')
        scoring_tasks = []

        path_list = os.listdir(path)

        for entry in path_list:
            full_path = os.path.join(path, entry)
            if os.path.isdir(full_path):
                if entry.startswith(ScoringTaskFactory.STAGE_FOUR_PREFIX):
                    if entry.endswith(ScoringTaskFactory.WEB_DATA_SUFFIX):
                        logger.debug('Skipping '+ entry + ' due to suffix')
                        continue

                    # we have a valid docking path
                    docktask = D3RTask(path, self.get_args())
                    docktask.set_stage(ScoringTaskFactory.DOCKSTAGE)
                    docktask.set_name(entry[
                                      len(ScoringTaskFactory.STAGE_FOUR_PREFIX)
                                      + 1:])
                    stask = ScoringTask(path,
                                        docktask.get_name() + '.' +
                                        ScoringTaskFactory.SCORING_SUFFIX,
                                        self.get_args())
                    if stask.can_run():
                        logger.debug('Adding task ' + stask.get_name())
                        scoring_tasks.append(stask)
                    else:
                        logger.debug(stask.get_name() + ' cannot be' +
                                     ' added : ' + stask.get_error())

        return scoring_tasks


class ScoringTask(D3RTask):
    """Performs Scoring

    """

    def __init__(self, path, name, docktask, args):
        super(ScoringTask, self).__init__(path, args)
        self.set_name(name)
        self.set_stage(5)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._docktask = docktask

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the docking task
           has `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies this task does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        self._docktask.update_status_from_filesystem()
        if self._docktask.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + self._docktask.get_name() + 'task' +
                        'has a status of ' + self._docktask.get_status())
            self.set_error(self._docktask.get_name() + ' task has ' +
                           self._docktask.get_status() + ' status')
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
        """Runs ScoringTask after verifying dock was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes scoring script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(ScoringTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('scoring set to ' +
                         self.get_args().scoring)
        except AttributeError:
            self.set_error('scoring not set')
            self.end()
            return

        try:
            logger.debug('pdbdb set to ' +
                         self.get_args().pdbdb)
        except AttributeError:
            self.set_error('pdbdb not set')
            self.end()
            return

        #
        # --pdbdb <path to pdb.extracted> --dockdir <stage.4.glide> \
        # --outdir <path to stage.5.glide.scoring>
        #
        cmd_to_run = (self.get_args().scoring + ' --pdbdb ' +
                      self.get_args().pdbdb + ' --dockdir ' +
                      self._docktask.get_dir() +
                      ' --outdir ' + self.get_dir())

        scoring_name = os.path.basename(self.get_args().scoring)

        self.run_external_command(scoring_name, cmd_to_run,
                                  True)
        # assess the result
        self.end()