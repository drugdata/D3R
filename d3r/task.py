# -*- coding: utf-8 -*-

import os
import logging
import re

logger = logging.getLogger(__name__)


class D3RParameters(object):
    """Holds parameters common to Tasks

    """
    pass


class UnsetPathException(Exception):
    """Exception to denote path is unset
    """
    pass


class UnsetNameException(Exception):
    """Exception to denote name is unset
    """
    pass


class UnsetStageException(Exception):
    """Exception to denote stage is unset
    """
    pass


class UnsetBlastDirException(Exception):
    """Exception to denote blastdir in D3RParameters is unset
    """
    pass


class TaskUnableToStartException(Exception):
    """Exception to denote when task cannot start due to failure
    """
    pass


class TaskFailedException(Exception):
    """Exception to denote when a task failed
    """
    pass


class D3RTask(object):
    """Represents a base Task that can be run.

    This is a base class from which other classes that actual do work
    can be derived from.  It provides internal attributes for name,
    stage, and status.  As well as a run() function.

    """
    STAGE_DIRNAME_PREFIX = "stage"

    START_FILE = "start"
    ERROR_FILE = "error"
    COMPLETE_FILE = "complete"

    START_STATUS = "start"
    COMPLETE_STATUS = "complete"
    UNKNOWN_STATUS = "unknown"
    NOTFOUND_STATUS = "notfound"
    ERROR_STATUS = "error"

    def __init__(self, path, args):
        """Constructor

        Creates a `D3RTask` with `D3RTask.UNKNOWN_STATUS` status
        with `path` and `args` set to values passed in
        """
        self._path = path
        self._name = None
        self._stage = None
        self._status = D3RTask.UNKNOWN_STATUS
        self._error = None
        self._args = args
        self._can_run = False

    def get_args(self):
        return self._args

    def set_args(self, args):
        self._args = args

    def get_path(self):
        return self._path

    def set_path(self, path):
        self._path = path

    def get_name(self):
        return self._name

    def set_name(self, name):
        self._name = name

    def get_stage(self):
        return self._stage

    def set_stage(self, stage):
        self._stage = stage

    def get_status(self):
        return self._status

    def set_status(self, status):
        self._status = status

    def get_error(self):
        return self._error

    def set_error(self, error):
        self._error = error

    def start(self):
        logger.info(self._name + ' task has started ')
        if not self.create_dir():
            self.set_error('Unable to create directory')
            self.end()
            return

        open(os.path.join(self._path, self.get_dir_name(),
                          D3RTask.START_FILE), 'a').close()

    def end(self):
        logger.info(self._name + ' task has finished with status ' +
                    self.get_status())
        if self._error is not None:
            logger.error(self._name + ' task failed with error ' +
                         self.get_error())
            open(os.path.join(self._path, self.get_dir_name(),
                              D3RTask.ERROR_FILE), 'a').close()

        if self.get_status() == D3RTask.COMPLETE_STATUS:
            open(os.path.join(self._path, self.get_dir_name(),
                              D3RTask.COMPLETE_FILE), 'a').close()

    def run(self):
        logger.info(self._name + ' task is running')

    def get_dir_name(self):
        """Gets directory name for task

           """
        if self._stage is None:
            raise UnsetStageException('Stage must be set')

        if self._name is None:
            raise UnsetNameException('Name must be set')

        return (D3RTask.STAGE_DIRNAME_PREFIX + "." + str(self._stage) +
                "." + self._name)

    def update_status_from_filesystem(self):
        """Updates status by querying filesystem.

           Sets status based on contents of path on filesystem.
           If path does not exist then status is NOTFOUND_STATUS.
           If complete file exists under path then status is COMPLETE_STATUS
           If error file exists under path then status is ERROR_STATUS
           If start file exists under path then status is START_STATUS
           else status is UNKNOWN_STATUS
           """
        if self._path is None:
            raise UnsetPathException('Path must be set')

        pathToCheck = os.path.join(self._path, self.get_dir_name())

        self.set_status(_get_status_of_task_in_dir(pathToCheck))
        logger.debug(self.get_name() + ' task status set to ' +
                     self.get_status())
        return self.get_status()

    def create_dir(self):
        """Creates directory for task.

           Directory will be named stage.<stage>.<name>
           and located under get_path()
           """
        if self.get_path() is None:
            logger.warning("Path is null cannot create directory")
            return

        if self.get_dir_name() is None:
            logger.warning("Dir name is null cannot create directory")
            return

        thePath = os.path.join(self.get_path(),
                               self.get_dir_name())

        logger.debug('Creating directory: ' + thePath)

        os.mkdir(thePath)

        if not os.path.isdir(thePath):
            logger.warning("Unable to create directory: " + thePath)
            return

        return thePath


class DataImportTask(D3RTask):
    """Represents DataImport Task

       """

    def __init__(self, path, args):
        super(DataImportTask, self).__init__(path, args)
        self.set_name('dataimport')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)


class MakeBlastDBTask(D3RTask):
    """Represents Blastable database

       This object does not follow the standard D3RTask
       structure.  This is because the blastable database
       is created externally in a legacy structure.

       """

    def __init__(self, path, args):
        """Constructor

           Constructor sets name to makeblastdb and ignores
           path set in <path> variable and instead sets
           path to args.blastdir. stage is set to 1

        """
        self.set_args(args)
        self.set_path(args.blastdir)
        self.set_name('makeblastdb')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def create_dir(self):
        """Creates path set in get_path()

           """
        os.makedirs(os.path.join(self._path, 'current'))

    def get_dir_name(self):
        """Will always return current

           """
        return 'current'


class BlastNFilterTask(D3RTask):
    """Performs Blast and filter of sequences

    """

    def __init__(self, path, args):
        super(BlastNFilterTask, self).__init__(path, args)
        self.set_name('blastnfilter')
        self.set_stage(2)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def can_run(self):
        self._can_run = False
        self._error = None
        # check blast
        makeblastdb = MakeBlastDBTask(self._path, self._args)
        makeblastdb.update_status_from_filesystem()
        if makeblastdb.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + makeblastdb.get_name() + 'task' +
                        'has a status of ' + makeblastdb.get_status())
            if makeblastdb.get_status() == D3RTask.ERROR_STATUS:
                self.set_error(makeblastdb.get_name() + ' task has ' +
                               makeblastdb.get_status() + ' status')
            return False

        # check data import
        dataImport = DataImportTask(self._path, self._args)
        dataImport.update_status_from_filesystem()
        if dataImport.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + dataImport.get_name() + 'task' +
                        'has a status of ' + dataImport.get_status())
            if dataImport.get_status() == D3RTask.ERROR_STATUS:
                self.set_error(dataImport.get_name() + ' task has ' +
                               dataImport.get_status() + ' status')
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
           creates a directory and invokes blastnfilter script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """

        if self._can_run is None:
            logger.info("Running can_run() to check if its allright " +
                        "to run")
            self.can_run()

        if self._can_run is False:
            logger.info("can_run() came back with false cannot run")
            return

        self.start()

        if self.get_error() is not None:
            self.end()
            return

        # Run the blastnfilter
        logger.info("Running blastnfilter")
        print "Running blastnfilter <> <>"

        self.set_status(D3RTask.COMPLETE_STATUS)
        # assess the result
        self.end()


def find_latest_year(celppdir):
    """Given a directory find the latest year

       The latest year will be a folder of 4 digits
       and have the highest value ie 2015
       :return: Directory or None if none is found
       :raises: Exception: If celppdir is not a directory
       """

    if not os.path.isdir(celppdir):
        raise Exception(celppdir+" is not a directory")

    dirPattern = re.compile("^\d\d\d\d$")

    latestYear = -1
    latestEntry = None
    fullPath = None
    for entry in os.listdir(celppdir):
        if re.match(dirPattern, entry):
            entryYear = int(entry)

            if entryYear > latestYear:
                fullPath = os.path.join(celppdir, entry)
                if os.path.isdir(fullPath):
                    latestYear = entryYear
                    latestEntry = fullPath

    if latestYear == -1:
        return

    return latestEntry


def find_latest_weekly_dataset(celppdir):
    """Given a directory find the latest dataset

       This method looks in directory passed in for
       paths with dataset.week.# format and returns
       the path with highest # on it.  If none are
       found then None is returned.
       :return: Directory upon success or None if no directory found
       :raises: Exception: If celppdir is not a directory
       """

    latestYear = find_latest_year(celppdir)

    if latestYear is None:
        return

    dirPattern = re.compile("^dataset.week.\d+$")

    latestEntry = None
    latestWeekNo = -1
    fullPath = None
    for entry in os.listdir(latestYear):
        if re.match(dirPattern, entry):
            weekNo = re.sub("dataset.week.", "", entry)
            if weekNo > latestWeekNo:
                fullPath = os.path.join(latestYear, entry)
                if os.path.isdir(fullPath):
                    latestWeekNo = weekNo
                    latestEntry = fullPath

    return latestEntry


def _get_status_of_task_in_dir(path):
    """Gets status of task based on existance of files in path

       Examines `path` and returns status based on the following
       conditions
       """
    if not os.path.isdir(path):
        return D3RTask.NOTFOUND_STATUS

    completeFile = os.path.join(path, D3RTask.COMPLETE_FILE)

    if os.path.isfile(completeFile):
        return D3RTask.COMPLETE_STATUS

    errorFile = os.path.join(path, D3RTask.ERROR_FILE)

    if os.path.isfile(errorFile):
        return D3RTask.ERROR_STATUS

    startFile = os.path.join(path, D3RTask.START_FILE)

    if os.path.isfile(startFile):
        return D3RTask.START_STATUS

    return D3RTask.UNKNOWN_STATUS
