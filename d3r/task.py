# -*- coding: utf-8 -*-

import os
import logging
import subprocess
import shlex
logger = logging.getLogger(__name__)


class D3RParameters(object):
    """Holds parameters common to Tasks

    """
    pass


class UnsetPathError(Exception):
    """Exception to denote path is unset
    """
    pass


class UnsetFileNameError(Exception):
    """Exception to denote file name is unset
    """
    pass


class UnsetNameError(Exception):
    """Exception to denote name is unset
    """
    pass


class UnsetStageError(Exception):
    """Exception to denote stage is unset
    """
    pass


class UnsetBlastDirError(Exception):
    """Exception to denote blastdir in D3RParameters is unset
    """
    pass


class TaskUnableToStartError(Exception):
    """Exception to denote when task cannot start due to failure
    """
    pass


class TaskFailedError(Exception):
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

        open(os.path.join(self.get_dir(),
                          D3RTask.START_FILE), 'a').close()

    def end(self):
        logger.info(self._name + ' task has finished with status ' +
                    self.get_status())
        if self._error is not None:
            logger.error(self._name + ' task failed with error ' +
                         self.get_error())
            open(os.path.join(self.get_dir(),
                              D3RTask.ERROR_FILE), 'a').close()

        if self.get_status() == D3RTask.COMPLETE_STATUS:
            open(os.path.join(self.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()

    def run(self):
        logger.info(self._name + ' task is running')

    def get_dir_name(self):
        """Gets directory name for task

           :raises: UnsetStageError, UnsetNameError
        """
        if self._stage is None:
            raise UnsetStageError('Stage must be set')

        if self._name is None:
            raise UnsetNameError('Name must be set')

        return (D3RTask.STAGE_DIRNAME_PREFIX + "." + str(self._stage) +
                "." + self._name)

    def get_dir(self):
        """Gets full path of Task

           :raises: UnsetPathError if path is not set,
                    and all exceptions from get_dir_name()
        """
        if self.get_path() is None:
            raise UnsetPathError('Path must be set')

        return os.path.join(self.get_path(),
                            self.get_dir_name())

    def update_status_from_filesystem(self):
        """Updates status by querying filesystem.

           Sets status based on contents of path on filesystem.
           If path does not exist then status is NOTFOUND_STATUS.
           If complete file exists under path then status is COMPLETE_STATUS
           If error file exists under path then status is ERROR_STATUS
           If start file exists under path then status is START_STATUS
           else status is UNKNOWN_STATUS
           """
        if self.get_path() is None:
            raise UnsetPathError('Path must be set')

        path_to_check = self.get_dir()

        self.set_status(self._get_status_of_task_in_dir(path_to_check))
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

        the_path = self.get_dir()

        logger.debug('Creating directory: ' + the_path)

        os.mkdir(the_path)

        if not os.path.isdir(the_path):
            logger.warning("Unable to create directory: " + the_path)
            return

        return the_path

    def write_to_file(self, str, file_name):
        """Writes `str` to `file_name` under `task` directory

           If `str` is None file is created but nothing is written.
           File is written using 'w' mode so if file exists it will
           be overwritten.

           :raises: UnsetFileNameError if file_name is None
        """
        if file_name is None:
            raise UnsetFileNameError('file_name must be set')
        file_to_write = os.path.join(self.get_dir(), file_name)
        logger.debug('Writing file ' + file_to_write)
        f = open(file_to_write, 'w')

        if str is not None:
            f.write(str)

        f.close()

    def _get_status_of_task_in_dir(self, path):
        """Gets status of task based on existance of files in path

           Examines get_path() and returns status based on the following
           conditions
        """
        if not os.path.isdir(path):
            return D3RTask.NOTFOUND_STATUS

        if os.path.isfile(os.path.join(path, D3RTask.COMPLETE_FILE)):
            return D3RTask.COMPLETE_STATUS

        if os.path.isfile(os.path.join(path, D3RTask.ERROR_FILE)):
            return D3RTask.ERROR_STATUS

        if os.path.isfile(os.path.join(path, D3RTask.START_FILE)):
            return D3RTask.START_STATUS

        return D3RTask.UNKNOWN_STATUS


class DataImportTask(D3RTask):
    """Represents DataImport Task

       """
    NONPOLYMER_TSV = "new_release_structure_nonpolymer.tsv"
    SEQUENCE_TSV = "new_release_structure_sequence.tsv"

    def __init__(self, path, args):
        super(DataImportTask, self).__init__(path, args)
        self.set_name('dataimport')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def get_nonpolymer_tsv(self):
        """Returns path to new_release_structure_nonpolymer.tsv file"""
        return os.path.join(self.get_dir() +
                            DataImportTask.NONPOLYMER_TSV)

    def get_sequence_tsv(self):
        """Returns path to new_release_structure_sequence.tsv file"""
        return os.path.join(self.get_dir() +
                            DataImportTask.SEQUENCE_TSV)


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
        """Determines if task can actually run

           This method first verifies the `MakeBlastDBTask` task
           and `DataImportTask` both have `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies a `BlastNFilterTask` does
           not already exist.
        """
        self._can_run = False
        self._error = None
        # check blast
        make_blastdb = MakeBlastDBTask(self._path, self._args)
        make_blastdb.update_status_from_filesystem()
        if make_blastdb.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + make_blastdb.get_name() + 'task' +
                        'has a status of ' + make_blastdb.get_status())
            if make_blastdb.get_status() == D3RTask.ERROR_STATUS:
                self.set_error(make_blastdb.get_name() + ' task has ' +
                               make_blastdb.get_status() + ' status')
            return False

        # check data import
        data_import = DataImportTask(self._path, self._args)
        data_import.update_status_from_filesystem()
        if data_import.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + data_import.get_name() + 'task' +
                        'has a status of ' + data_import.get_status())
            if data_import.get_status() == D3RTask.ERROR_STATUS:
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

        data_import = DataImportTask(self._path, self._args)

        cmd_to_run = (self.get_args().blastnfilter + ' --nonpolymertsv ' +
                      data_import.get_nonpolymer_tsv() +
                      ' --sequencetsv ' +
                      data_import.get_sequence_tsv() +
                      ' --outdir ' + self.get_dir())

        # Run the blastnfilter
        logger.info("Running command " + cmd_to_run)
        try:
            p = subprocess.Popen(shlex.split(cmd_to_run),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        except Exception as e:
            logger.exception("Error caught exception")
            self.set_status(D3RTask.ERROR_STATUS)
            self.set_error("Caught Exception trying to run " +
                           cmd_to_run + " : " + e.message)
            self.end()
            return

        out, err = p.communicate()

        self.write_to_file(err,
                           os.path.basename(self.get_args().blastnfilter +
                                            '.stderr'))
        self.write_to_file(out,
                           os.path.basename(self.get_args().blastnfilter +
                                            '.stdout'))

        if p.returncode == 0:
            self.set_status(D3RTask.COMPLETE_STATUS)
        else:
            self.set_status(D3RTask.ERROR_STATUS)
            self.set_error("Non zero exit code: " + p.returncode +
                           "received. Standard out: " + out +
                           " Standard error : " + err)

        # assess the result
        self.end()
