# -*- coding: utf-8 -*-

import os
import time
import logging
import subprocess
import shlex
import smtplib
import platform
from email.mime.text import MIMEText

from d3r.celpp.uploader import FtpFileUploader


logger = logging.getLogger(__name__)


class D3RParameters(object):
    """Holds parameters common to Tasks

    """
    pass


class TaskException(Exception):
    """base exception class for all the other exceptions provided by
       this module.
    """
    pass


class UnsetPathError(TaskException):
    """Exception to denote path is unset
    """
    pass


class UnsetFileNameError(TaskException):
    """Exception to denote file name is unset
    """
    pass


class UnsetNameError(TaskException):
    """Exception to denote name is unset
    """
    pass


class UnsetCommandError(TaskException):
    """Exception to denote command is unset
    """
    pass


class UnsetStageError(TaskException):
    """Exception to denote stage is unset
    """
    pass


class TaskUnableToStartError(TaskException):
    """Exception to denote when task cannot start due to failure
    """
    pass


class TaskFailedError(TaskException):
    """Exception to denote when a task failed
    """
    pass


class D3RTask(object):
    """Represents a base Task that can be run.

    This is a base class from which other classes that actual do work
    can be derived from.  It provides internal attributes for name,
    stage, and status.  As well as a run() function.

    """

    # Put at start of subject line in all emails sent from D3RTask
    SUBJECT_LINE_PREFIX = "[d3rcelpp] "

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
        self._can_run = None
        self._start_time = None
        self._duration = -1
        self._email_log = None

        try:
            logger.debug('ftpconfig set to ' + args.ftpconfig)
            self._file_uploader = FtpFileUploader(args.ftpconfig)
        except AttributeError:
            logger.debug('ftpconfig not set')
            self._file_uploader = None

    def get_args(self):
        return self._args

    def set_args(self, args):
        self._args = args

    def set_email_log(self, email_log):
        self._email_log = email_log

    def append_to_email_log(self, append_log):
        if self._email_log is None:
            self._email_log = append_log
        else:
            self._email_log += append_log

    def get_email_log(self):
        return self._email_log

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

    def get_file_uploader(self):
        """Gets the file uploader
        """
        return self._file_uploader

    def set_file_uploader(self, file_uploader):
        """Sets file uploader
        """
        self._file_uploader = file_uploader

    def get_uploadable_files(self):
        """Returns a list of files that can be uploaded to remote server

           Provides a way for a D3RTask to define what files should
           uploaded to remote server by any D3RTaskUploader.  The default
           implementation returns an empty list
           :returns: list of files that can be uploaded.  The files should
           have full paths
        """
        return []

    def _get_smtp_server(self):
        """Gets smtplib server for localhost

           creates smtplib.SMTP server using get_args().smtp and
           get_args().smtpport for the smtp server and port
        :return: SMTP server
        """
        return smtplib.SMTP(self.get_args().smtp,
                            self.get_args().smtpport)

    def _build_from_address(self):
        """Returns from email address

           :return: string in format login@host
        """
        hostname = platform.node()
        if hostname is '':
            hostname = 'localhost'
        return (os.getlogin() + '@' + hostname)

    def _get_program_name(self):
        """Gets name of program running

        Looks at _args.program for the program name, if unset or
        empty then __file__ is used
        :return: Name of program that invoked this
        """
        program_name = __file__
        version = ''
        try:
            program_name = self.get_args().program
            version = self.get_args().version
        except:
            logger.debug('args.program or args.version was not set.  '
                         'using __file__')
        return program_name + ' ' + version

    def _get_time(self):
        """ Gets the current time as a string using ctime()
        :return: ctime string followed by --
        """
        return time.ctime() + ' -- '

    def _send_start_email(self):
        """Creates start email message that can be passed to sendmail

        """
        email_log = ''
        if self._email_log is not None:
            email_log = (self._email_log + '\n\n')

        the_message = ('Hi,\n ' + self.get_dir_name() + ' has started' +
                       ' running on host ' + platform.node() +
                       '\nunder path ' + self.get_dir() + '\n\n' + email_log +
                       'Sincerely,\n\n' + self._get_program_name())

        msg = MIMEText(the_message)
        msg['Subject'] = (D3RTask.SUBJECT_LINE_PREFIX + self._get_time() +
                          self.get_dir_name() +
                          ' has started running')
        self._send_email(msg)

    def _send_end_email(self):
        """Creates end email message that can be passed to sendmail

        """

        email_log = ''
        if self._email_log is not None:
            email_log = (self._email_log + '\n\n')

        error_msg = ''
        if self.get_error() is not None:
            error_msg = ('Error:\n' + self.get_error() + '\n\n')

        the_message = ('Hi,\n ' + self.get_dir_name() + ' has finished' +
                       ' after ' + str(self._duration) + ' seconds' +
                       ' with status ' + self.get_status() + ' on host '
                       + platform.node() +
                       '\nunder path ' + self.get_dir() + '\n\n' + email_log +
                       error_msg + 'Sincerely,\n\n' + self._get_program_name())

        msg = MIMEText(the_message)
        msg['Subject'] = (D3RTask.SUBJECT_LINE_PREFIX + self._get_time() +
                          self.get_dir_name() +
                          ' has finished with status ' + self.get_status())
        self._send_email(msg)

    def _send_email(self, mime_msg):
        """Sends email with message passed in


           If email attribute is set with valid email addresses
           in get_args() object then this method will send
           an email with `message` passed in as content.
        """
        try:
            if self.get_args().email is not None:
                email_list = self.get_args().email.split(',')
                logger.debug('Sending email to: ' + ", ".join(email_list))
                server = self._get_smtp_server()

                from_addr = self._build_from_address()

                mime_msg['From'] = from_addr
                mime_msg['To'] = ", ".join(email_list)

                if self.get_args().loglevel == 'DEBUG':
                    server.set_debuglevel(1)

                server.sendmail(from_addr,
                                email_list, mime_msg.as_string())
                server.quit()
                return
            else:
                logger.debug('Email was set to None, skipping email ' +
                             'notification')
        except AttributeError as e:
            logger.debug('No email set, skipping email notification ' +
                         str(e))

    def _upload_task(self):
        """Uploads files related to task with FileUploader passed in


           :return: None
        """
        if self.get_file_uploader() is None:
            return
        try:
            uploadable_files = self.get_uploadable_files()
            self.get_file_uploader().upload_files(uploadable_files)
            summary = self.get_file_uploader().get_upload_summary()
            self.append_to_email_log('\n' + summary + '\n')
        except:
            logger.exception('Caught exception trying to upload files')


    def start(self):
        """Denotes start of task

           Creates task directory, along with token start file.
           If get_args() has email attribute then an email is
           sent denoting task start.  Also sets get_status()
           to D3RTask.START_STATUS.
           If there is an error creating task directory then
           end() is invoked and get_error() will be set with
           message denoting error
        """
        self._start_time = int(time.time())
        logger.info(self.get_dir_name() + ' task has started ')

        self.set_status(D3RTask.START_STATUS)

        self._send_start_email()

        try:
            self.create_dir()
        except Exception as e:
            self.set_error('Caught exception trying to create directory ' +
                           str(e))
            self.end()
            return

        open(os.path.join(self.get_dir(),
                          D3RTask.START_FILE), 'a').close()

    def end(self):
        """Denotes task has completed

           If get_error() is not none then 'error' token file is
           put in task directory and get_status() set to D3RTask.ERROR_STATUS
           Otherwise 'complete' token file is put in task directory.
           if get_args() email is set then a notification of task
           completion is sent
        """

        if self._start_time is not None:
            self._duration = int(time.time()) - self._start_time
        else:
            self._duration = -1

        logger.info(self.get_dir_name() + ' task took ' +
                    str(self._duration) + ' seconds and has ' +
                    'finished with status ' + self.get_status())

        if self.get_error() is not None:
            logger.error(self.get_dir() + ' task failed with error ' +
                         self.get_error())
            self.set_status(D3RTask.ERROR_STATUS)
            try:
                open(os.path.join(self.get_dir(),
                                  D3RTask.ERROR_FILE), 'a').close()
            except:
                logger.exception('Caught exception')

        self._upload_task()

        if self.get_status() == D3RTask.ERROR_STATUS:
            try:
                open(os.path.join(self.get_dir(),
                                  D3RTask.ERROR_FILE), 'a').close()
            except:
                logger.exception('Caught Exception')
            self._send_end_email()
            return

        self.set_status(D3RTask.COMPLETE_STATUS)
        open(os.path.join(self.get_dir(),
                          D3RTask.COMPLETE_FILE), 'a').close()
        self._send_end_email()

    def can_run(self):
        """Always returns False

           Derived classes should override this method
        """
        self._can_run = False
        return False

    def run(self):
        """Sees if task can run

           This is the base implementation and does not actually
           run anything.  Instead this method calls can_run()
           if the internal variable _can_run is None.  It is assumed
           can_run() will set _can_run variable.  If _can_run is False
           method will just return.  If _can_run is True then start()
           is called and after invocation get_error() is checked, if
           an error is found end() is called and _can_run is set to False
        """
        name = self._name
        if name is None:
            name = 'unset'

        logger.info(name + ' task is running')
        if self._can_run is None:
            logger.info("Running can_run() to check if its allright " +
                        "to run")
            self.can_run()

        if self._can_run is False:
            logger.info("can_run() came back with false cannot run")
            if self.get_error() is not None:
                self.end()
            return

        self.start()

        if self.get_error() is not None:
            self._can_run = False
            self.end()
            return

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
           :return: status as string updated from file system
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
           :raises: UnsetPathError if path is not set and all exceptions
                    from get_dir_name(), and OSError if there is a
                    problem creating the directory
           :return: path to created directory
           """
        the_path = self.get_dir()

        logger.debug('Creating directory: ' + the_path)

        os.mkdir(the_path)

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
           :return: status of task as string
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

    def run_external_command(self, command_name, cmd_to_run,
                             command_failure_is_fatal):
        """Runs external command line process

        Method runs external process logging the command
        run to email log.  Standard output and error are
        written to current working directory `command_name`.stdout
        and `command_name`.stderr respectively and if process.

        :param command_name: Name of command, human readable version
                             used as prefix for .stderr and .stdout
                             files
        :param cmd_to_run: command with arguments to run,
                        passed to subprocess.Popen
        :param command_failure_is_fatal: True means if process
               fails with nonzero exit code or there is an
               exception then status of task is set to
               D3RTask.ERROR_STATUS and failure message appended
               to email log and task.set_error is set
        :return: exit code of process
        :raise: UnsetNameError if command_name is None
        :raise: UnsetCommandError if cmd_to_run is None
        """
        if command_name is None:
            raise UnsetNameError('Command name must be set')

        if cmd_to_run is None:
            raise UnsetCommandError('Command must be set')

        if command_failure_is_fatal is None:
            command_failure_is_fatal = True

        logger.info("Running command " + cmd_to_run)

        self.append_to_email_log('Running command: ' + cmd_to_run + '\n')
        try:
            p = subprocess.Popen(shlex.split(cmd_to_run),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        except Exception as e:
            logger.exception("Error caught exception")
            self.set_status(D3RTask.ERROR_STATUS)
            self.set_error("Caught Exception trying to run " +
                           cmd_to_run + " : " + str(e))
            self.end()
            return 1

        out, err = p.communicate()

        self.write_to_file(err, command_name + '.stderr')
        self.write_to_file(out, command_name + '.stdout')

        if p.returncode != 0:
            if command_failure_is_fatal:
                self.set_status(D3RTask.ERROR_STATUS)
                self.set_error("Non zero exit code: " + str(p.returncode) +
                               " received. Standard out: " + out +
                               " Standard error: " + err)
            else:
                self.append_to_email_log("Although considered non fatal " +
                                         "for processing of stage a " +
                                         "non zero exit code: " +
                                         str(p.returncode) +
                                         "received. Standard out: " + out +
                                         " Standard error : " + err)
        return p.returncode
