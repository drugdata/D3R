# -*- coding: utf-8 -*-

import os
import pwd
import time
import logging
import smtplib
import platform
import mimetypes
import shutil
import uuid
import configparser
from email import encoders
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from d3r.celpp.filetransfer import FtpFileTransfer
from d3r.celpp import util


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


class EmailSendError(TaskException):
    """Exception to denote when SmtpEmailer was unable to send an email
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
    STDERR_SUFFIX = ".stderr"
    STDOUT_SUFFIX = ".stdout"
    TMP_DIR_SUFFIX = '.tmpdir'
    MAX_CHARS_FOR_EMAIL_STR = 250000
    TEXT_TRUNCATED_STR = 'TEXT TRUNCATED\n'

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
            self._file_uploader = FtpFileTransfer(args.ftpconfig)
        except Exception:
            logger.debug('FtpFileUploader not set.  This may not be '
                         'an error')
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

    def get_file_transfer(self):
        """Gets the file uploader
        """
        return self._file_uploader

    def set_file_transfer(self, file_uploader):
        """Sets file uploader
        """
        self._file_uploader = file_uploader

    def get_version_from_start_file(self):
        """Gets version of CELPP/D3R run by parsing start file
           :returns version from start file or unknown if file is empty or
                    does not exist
        """
        version = 'unknown'
        sfile = os.path.join(self.get_dir(), self.START_FILE)
        if os.path.isfile(sfile):
            with open(sfile, 'r') as fp:
                data = fp.read()
                if len(data) > 0:
                    version = str(data)
        return version

    def get_uploadable_files(self):
        """Returns a list of files that can be uploaded to remote server

           Provides a way for a D3RTask to define what files should
           uploaded to remote server by any D3RTaskUploader.  The default
           implementation returns an empty list
           :returns: list of files that can be uploaded.  The files should
           have full paths
        """
        file_list = []

        try:
            out_dir = self.get_dir()
        except Exception:
            logger.exception('Unable to get dir for task')
            return file_list
        try:
            for entry in os.listdir(out_dir):
                if not os.path.isfile(os.path.join(out_dir, entry)):
                    continue

                if entry.endswith(D3RTask.STDERR_SUFFIX):
                    file_list.append(os.path.join(out_dir, entry))
                elif entry.endswith(D3RTask.STDOUT_SUFFIX):
                    file_list.append(os.path.join(out_dir, entry))

            error_file = os.path.join(out_dir, D3RTask.ERROR_FILE)
            if os.path.isfile(error_file):
                logger.debug('Found error file adding to list')
                file_list.append(error_file)

        except OSError:
            logger.warning('Caught exception looking for error and '
                           'stderr/stdout files')

        return file_list

    def _get_program_name(self):
        """Gets name of program running

        Looks at _args.program for the program name, if unset or
        empty then __file__ is used
        :return: Name of program that invoked this
        """
        program_name = __file__
        try:
            program_name = self.get_args().program
        except Exception:
            logger.debug('args.program  was not set.  '
                         'using __file__')
        return program_name + ' ' + self.get_program_version()

    def get_program_version(self):
        """Gets version of program running
        This is done by looking at _args.version
        :return: Version of program or empty string if unset
        """
        try:
            version = self.get_args().version
        except AttributeError:
            version = ''
            logger.debug('args.version was not set')
        return str(version)

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

        subject = (D3RTask.SUBJECT_LINE_PREFIX + self._get_time() +
                   self.get_dir_name() + ' has started running')
        try:
            self._send_email(self.get_args().email, subject, the_message)
        except AttributeError as e:
            logger.debug(str(e) + ' email not set. skipping')

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

        subject = (D3RTask.SUBJECT_LINE_PREFIX + self._get_time() +
                   self.get_dir_name() + ' has finished with status ' +
                   self.get_status())
        try:
            try:
                self._send_email(self.get_args().email, subject, the_message)
            except AttributeError as e:
                logger.debug(str(e) + ' email not set. skipping')

            if self.get_error() is not None:
                logger.debug("Sending summary email cause task failed")
                try:
                    self._send_email(self.get_args().summaryemail, subject,
                                     the_message)
                except AttributeError as e:
                    logger.debug(str(e) + ' summary email not set. skipping')

        except Exception:
            logger.exception('Caught exception trying to send end email, '
                             'skipping email notification')

    def _send_email(self, to_email_csv, subject, message):
        """Sends email with message passed in


           If email attribute is set with valid email addresses
           in get_args() object then this method will send
           an email with `message` passed in as content.
        """
        try:
            if to_email_csv is not None:
                fac = SmtpEmailerFactory(self.get_args())
                emailer = fac.get_smtp_emailer()

                split_email_list = to_email_csv.split(',')

                logger.debug('Sending email to: ' +
                             ", ".join(split_email_list))

                emailer.send_email(split_email_list, subject, message)
                return
            else:
                logger.debug('Email was set to None, skipping email ' +
                             'notification')
        except AttributeError as e:
            logger.debug('No email passed in, skipping email notification ' +
                         str(e))

    def _upload_task(self):
        """Uploads files related to task with FileUploader passed in


           :return: None
        """
        if self.get_file_transfer() is None:
            return
        try:
            uploadable_files = self.get_uploadable_files()
            self.get_file_transfer().connect()
            self.get_file_transfer().upload_files(uploadable_files)
            summary = self.get_file_transfer().get_upload_summary()
            self.append_to_email_log('\n' + summary + '\n')
        except Exception:
            logger.exception('Caught exception trying to upload files')
        finally:
            self.get_file_transfer().disconnect()

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

        sfile = open(os.path.join(self.get_dir(),
                                  D3RTask.START_FILE), 'w')
        sfile.write(self.get_program_version())
        sfile.flush()
        sfile.close()

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
                f = open(os.path.join(self.get_dir(),
                                      D3RTask.ERROR_FILE), 'a')
                f.write(self.get_error() + '\n')
                f.flush()
                f.close()
            except Exception:
                logger.exception('Caught exception')

        self._upload_task()

        if self.get_status() == D3RTask.ERROR_STATUS:
            try:
                err_file = os.path.join(self.get_dir(),
                                        D3RTask.ERROR_FILE)
                if not os.path.isfile(err_file):
                    open(err_file, 'a').close()
            except Exception:
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
        try:
            logger.debug('Checking path ' + path_to_check)
        except Exception:
            pass

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
                             command_failure_is_fatal,
                             timeout=None,
                             kill_delay=10,
                             polling_sleep_time=1):
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
        cmd_tmp_dir = None
        try:
            if timeout is not None:
                cmd_tmp_dir = os.path.join(self.get_dir(),
                                           command_name + str(uuid.uuid1()) +
                                           D3RTask.TMP_DIR_SUFFIX)
                os.makedirs(cmd_tmp_dir, mode=0o0775)
                pst = polling_sleep_time  # trying to make flake8 happy

                returncode, out, err = util.\
                    run_external_command_with_timeout(cmd_to_run,
                                                      cmd_tmp_dir,
                                                      timeout=timeout,
                                                      kill_delay=kill_delay,
                                                      polling_sleep_time=pst)
            else:
                returncode, out, err = util.run_external_command(cmd_to_run)
        except Exception as e:
            logger.exception("Error caught exception")
            self.set_status(D3RTask.ERROR_STATUS)
            self.set_error("Caught Exception trying to run " +
                           cmd_to_run + " : " + str(e))
            self.end()
            return 1
        finally:
            if cmd_tmp_dir is not None:
                shutil.rmtree(cmd_tmp_dir)

        self.write_to_file(err, command_name + D3RTask.STDERR_SUFFIX)
        self.write_to_file(out, command_name + D3RTask.STDOUT_SUFFIX)

        if returncode != 0:
            if command_failure_is_fatal:
                self.set_status(D3RTask.ERROR_STATUS)
                self.set_error("Non zero exit code: " + str(returncode) +
                               " received. Standard out: " + out +
                               " Standard error: " + err)
            else:
                max_chars = D3RTask.MAX_CHARS_FOR_EMAIL_STR
                trunc_out = self._get_email_truncated_string(out,
                                                             max_chars)
                trunc_err = self._get_email_truncated_string(err,
                                                             max_chars)

                self.append_to_email_log("Although considered non fatal " +
                                         "for processing of stage a " +
                                         "non zero exit code: " +
                                         str(returncode) +
                                         "received. Standard out: " +
                                         trunc_out +
                                         " Standard error : " + trunc_err)
        return returncode

    def _get_email_truncated_string(self, val, max_chars):
        """Truncates string so it fits within a certain number of characters
        This method will remove characters from the start of the string so
        it fits within the `max_chars` limit and prefixed with value of
        `D3RTask.TEXT_TRUNCATED_STR` to the start.
        NOTE: the prefix string is not counted toward the max chars
        :param val: String to truncate
        :params max_chars: Maximum length of string
        :return: Truncated `val` string with `D3RTask.TEXT_TRUNCATED_STR`
                 prefixed
        """
        if val is None:
            logger.warning('String passed in was None so returning None')
            return None

        if max_chars is None or max_chars < 0:
            logger.debug('max_chars is set to None or is negative, '
                         'just returning string')
            return val

        str_len = len(val)
        logger.debug('String length is: ' + str(str_len))
        if str_len <= max_chars:
            return val

        trunc_start = str_len - max_chars
        logger.debug('Truncating string at index: ' + str(trunc_start))
        return (D3RTask.TEXT_TRUNCATED_STR +
                val[trunc_start:])


class Attachment(object):
    """Class to hold file to attach to email
    """
    def __init__(self, file_to_attach, desired_name):
        self._file_to_attach = file_to_attach
        self._desired_name = desired_name

    def get_file_to_attach(self):
        return self._file_to_attach

    def get_desired_name(self):
        return self._desired_name


class BaseConfig(object):
    """Base class for ConfigParser objects containing functions
       usable by derived classes
    """
    def __init__(self):
        """Constructor"""
        pass

    def _get_config(self, configfile):
        if not os.path.isfile(configfile):
            logger.warning(configfile + ' is not a file')
            return None

        config = configparser.ConfigParser()
        config.read(configfile)
        return config

    def _get_value(self, config, section, option):
        """Calls get() on configparser object `config` passed in to
           get value for corresponding option
           :param config: ConfigParser object loaded with smtp configuration
           :param section: Section to look for value
           :param option: Keyword value or value left of =
                           (ie foo = X the keyword would be foo)
           :returns: value as string or whatever ConfigParser().get()
                     returns or None if not found
        """
        if section is None:
            logger.error('Section cannot be None')
            return None

        if not config.has_section(section):
            logger.warning(section +
                           ' section not found in configuration')
            return None

        if not config.has_option(section, option):
            logger.info('In parsing ' + section + ' configuration ' + option +
                        ' not found')
            return None

        return config.get(section, option)


class WebsiteServiceConfig(BaseConfig):
    """Parses configfile passed into constructor to
       obtain parameters to connect to website REST services
       The configfile should be in configparser format as seen
       in this example:

       [websiteservice]

       url = https://blah.com/api/1/d3r/celpp
       apikey = asdlkjas
       basicauth_user = joe
       basicauth_password = secret123
       source = dev
       portal_name = d3r
       timeout = 5
    """
    DEFAULT = 'websiteservice'
    WEB_URL = 'url'
    WEB_APIKEY = 'apikey'
    WEB_SOURCE = 'source'
    WEB_PORTAL_NAME = 'portal_name'
    WEB_BASIC_USER = 'basicauth_user'
    WEB_BASIC_PASS = 'basicauth_password'
    WEB_TIMEOUT = 'timeout'
    RMSD_ENDPOINT = 'rmsd'
    TARGETS_ENDPOINT = 'week'
    JSON_KEY = 'json'
    CONTENT_TYPE_VAL = 'application/json'
    CONTENT_TYPE_KEY = 'Content-Type'
    API_KEY_KEY = 'X-API-KEY'

    def __init__(self, configfile=None):
        """Parses configuration file set in `configfile`
        """
        super(BaseConfig, self).__init__()
        self._url = None
        self._apikey = None
        self._source = 'notset'
        self._portal_name = 'notset'
        self._basicauth_user = None
        self._basicauth_pass = None
        self._timeout = 0.1

        if configfile is not None:
            try:
                self._parse_config(configfile)
            except configparser.Error:
                logger.exception('Caught Exception trying to parse ' +
                                 configfile)

    def _parse_config(self, configfile):
        """Parses configuration file passed in
           to set internal variables such as user, host, password, etc.
        """
        config = self._get_config(configfile)
        if config is None:
            return

        self._url = self._get_value(config, WebsiteServiceConfig.DEFAULT,
                                    WebsiteServiceConfig.WEB_URL)
        self._apikey = self._get_value(config, WebsiteServiceConfig.DEFAULT,
                                       WebsiteServiceConfig.WEB_APIKEY)
        self._source = self._get_value(config, WebsiteServiceConfig.DEFAULT,
                                       WebsiteServiceConfig.WEB_SOURCE)
        self._portal_name = self._get_value(config,
                                            WebsiteServiceConfig.DEFAULT,
                                            WebsiteServiceConfig.
                                            WEB_PORTAL_NAME)

        self._basicauth_user = self._get_value(config,
                                               WebsiteServiceConfig.DEFAULT,
                                               WebsiteServiceConfig.
                                               WEB_BASIC_USER)
        self._basicauth_pass = self._get_value(config,
                                               WebsiteServiceConfig.DEFAULT,
                                               WebsiteServiceConfig.
                                               WEB_BASIC_PASS)

        timeout = self._get_value(config, WebsiteServiceConfig.DEFAULT,
                                  WebsiteServiceConfig.WEB_TIMEOUT)
        try:
            if timeout is not None:
                self._timeout = float(timeout)
        except SyntaxError:
            logger.error('Unable to convert timeout : ' + str(timeout) +
                         ' to float')
        except ValueError:
            logger.error('Unable to convert timeout : ' + str(timeout) +
                         ' to float')

    def get_rmsd_url(self):
        """Gets website REST service base url with rmsd endpoint
        """
        if self._url is None:
            return None
        if self._url.endswith('/'):
            return self._url + WebsiteServiceConfig.RMSD_ENDPOINT
        return self._url + '/' + WebsiteServiceConfig.RMSD_ENDPOINT

    def get_apikey(self):
        """Get api key
        """
        return self._apikey

    def get_targets_url(self):
        """Gets website REST service base url with targets aka week endpoint
        """
        if self._url is None:
            return None
        if self._url.endswith('/'):
            return self._url + WebsiteServiceConfig.TARGETS_ENDPOINT
        return self._url + '/' + WebsiteServiceConfig.TARGETS_ENDPOINT

    def get_source(self):
        """Gets source
        """
        return self._source

    def get_portal_name(self):
        """Gets portal_name
        """
        return self._portal_name

    def get_basicauth_user(self):
        """Gets basic auth username
        """
        return self._basicauth_user

    def get_basicauth_password(self):
        """Gets basic auth password
        """
        return self._basicauth_pass

    def get_timeout(self):
        """Gets web request timeout
        """
        return self._timeout


class SmtpConfig(BaseConfig):
    """Parses configfile passed into constructor to
       obtain parameters to connect to smtp server.
       The configfile should be in configparser format
       as seen in this example:

       [smtp]
       host = foo.com
       user = bob
       port = 25
       password = 12345
       from_address = no-reply@bob.com
       replyto_address = uh@bob.com
    """
    DEFAULT = 'smtp'
    SMTP_HOST = 'host'
    SMTP_PORT = 'port'
    SMTP_USER = 'user'
    SMTP_PASS = 'password'
    SMTP_FROM_ADDRESS = 'from_address'
    SMTP_REPLYTO_ADDRESS = 'replyto_address'
    DEFAULT_PORT = 25

    def __init__(self, configfile=None):
        """If set parses configuration file set in `configfile`
           setting smtp connection information. If
           port is not set in configuration, 25 is used  and
           if host is not set then localhost is used.

        :param configfile: path to config file containing
                           information to connect to smtp
                           server. See constructor
                           comments for structure of
                           configuration file

        :return: SmtpConfig object
        """
        super(BaseConfig, self).__init__()
        self._user = None
        self._host = 'localhost'
        self._port = SmtpConfig.DEFAULT_PORT
        self._password = None
        self._from_address = None
        self._replyto_address = None

        if configfile is not None:
            try:
                self._parse_config(configfile)
            except configparser.Error:
                logger.exception('Caught Exception trying to parse ' +
                                 configfile)

    def _parse_config(self, configfile):
        """Parses configuration file passed in
           to set internal variables such as user, host, password, etc.
        """
        config = self._get_config(configfile)
        if config is None:
            return

        self._user = self._get_value(config, SmtpConfig.DEFAULT,
                                     SmtpConfig.SMTP_USER)

        host = self._get_value(config, SmtpConfig.DEFAULT,
                               SmtpConfig.SMTP_HOST)
        if host is not None:
            self._host = host

        port = self._get_value(config, SmtpConfig.DEFAULT,
                               SmtpConfig.SMTP_PORT)
        if port is not None:
            try:
                self._port = int(port)
            except ValueError:
                logger.exception('Unable to convert port to int. '
                                 'Using default')
                self._port = SmtpConfig.DEFAULT_PORT

        self._password = self._get_value(config, SmtpConfig.DEFAULT,
                                         SmtpConfig.SMTP_PASS)

        self._from_address = self._get_value(config, SmtpConfig.DEFAULT,
                                             SmtpConfig.SMTP_FROM_ADDRESS)

        self._replyto_address = self._get_value(config, SmtpConfig.DEFAULT,
                                                SmtpConfig.
                                                SMTP_REPLYTO_ADDRESS)

    def get_host(self):
        """Gets host"""
        return self._host

    def get_port(self):
        """Gets port as int"""
        return self._port

    def get_user(self):
        """Gets user"""
        return self._user

    def get_password(self):
        """Gets password"""
        return self._password

    def get_from_address(self):
        """Gets from address"""
        return self._from_address

    def get_replyto_address(self):
        """Gets replyto address"""
        return self._replyto_address


class SmtpEmailerFactory(object):
    """Factory that parses command line arguments
       passed in via `args` in constructor to get
       connection information and create a SmtpEmailer
       object.
    """

    def __init__(self, args):
        """Constructor
        :param args: arguments set via argparse. This
                     constructor looks for args.smtpconfig
                     and if found it is assumed to be a path
                     to smtp configuration file. See SmtpConfig
                     class for the format of this configuration
                     file.
        """
        configfile = None
        if args is not None:
            try:
                configfile = args.smtpconfig
            except AttributeError:
                logger.info('No smtpconfig configuration file set')
        self._smtpconfig = SmtpConfig(configfile=configfile)

    def get_smtp_emailer(self):
        """Gets SmtpEmailer object
        """
        return SmtpEmailer(smtp_host=self._smtpconfig.get_host(),
                           port=self._smtpconfig.get_port(),
                           user=self._smtpconfig.get_user(),
                           password=self._smtpconfig.get_password(),
                           fromaddr=self._smtpconfig.get_from_address(),
                           replyto=self._smtpconfig.get_replyto_address())


class SmtpEmailer(object):
    """Simple Wrapper class to send email via smtplib
    """

    def __init__(self, smtp_host='localhost', port=25, user=None,
                 password=None, fromaddr=None, replyto=None):
        """Constructor
        :param smtp_host: host of stmp server
        :param port: port stmp server is on
        """
        self._smtp_host = smtp_host
        self._port = port
        self._alt_smtp_server = None
        self._user = user
        self._password = password
        if fromaddr is None:
            self._fromaddr = self._generate_from_address_using_login_and_host()
        else:
            self._fromaddr = fromaddr

        self._replyto = replyto

    def set_alternate_smtp_server(self, server):
        """Sets alternate smtp server to use
        """
        self._alt_smtp_server = server

    def send_email(self, to_list, subject, message, attachments=None):
        """Sends email
        :param from_address: from email address
        :param to_list: list of email addresses as strings to send email to
        :param subject: Subject of email
        :param message: Body of email
        :param reply_to: optional parameter to alter the reply to email
        :param attachments: optional list parameter containing Attachment
        objects
        representing files to attach to email
        :param htmlmessage: optional parameter containing html version of
        email body
        :raises EmailSendError: if there was an error creating or sending
        the email
        """
        server = None
        try:

            mime_msg = self._build_mime_message(self._fromaddr, to_list,
                                                subject, message,
                                                self._replyto, attachments)
            server = self._get_server()
            server.sendmail(self._fromaddr, to_list, mime_msg.as_string())
        except smtplib.SMTPConnectError as e:
            logger.exception('Caught exception')
            raise EmailSendError('Unable to connect to smtp server ' + str(e))
        except Exception as e:
            logger.exception('Caught exception')
            raise EmailSendError('Caught exception ' + str(e))
        finally:
            if server is not None:
                server.quit()

    def _generate_from_address_using_login_and_host(self, hostname=None):
        """Creates from email address from login and hostname
        of machine running this script.
        :returns: from email address as string
        """
        if hostname is None:
            hostname = platform.node()

        if hostname is '':
            hostname = 'localhost'
        return pwd.getpwuid(os.getuid())[0] + '@' + hostname

    def _get_server(self):
        """Gets Smtp server
        :returns: smtplib.SMTP object unless set_alternate_smtp_server
        was set in
        which case that object is returned
        """
        if self._alt_smtp_server is not None:
            smtp = self._alt_smtp_server
        else:
            smtp = smtplib.SMTP(self._smtp_host,
                                self._port)
        if self._password is not None and self._user is not None:
            logger.debug('Attempting SMTP login with user: ' + str(self._user))
            smtp.login(self._user, self._password)
        return smtp

    def _append_attachments(self, msg_root, attachments):
        """Attaches attachments to `msg_root`
        :param msg_root: MIMEMultipart message to attach the files to
        :param attachments: List of Attachment objects to attach to msg_root
        :returns msg_root: with attachments added
        """
        if attachments is None:
            return msg_root

        for a in attachments:
            fname = a.get_file_to_attach()
            if not os.path.isfile(fname):
                logger.error('File does not exist: ' + fname)
                continue
            ctype, encoding = mimetypes.guess_type(fname)
            if ctype is None or encoding is not None:
                # go with binary type if we did not get a guess
                ctype = 'application/octet-stream'
            maintype, subtype = ctype.split('/', 1)
            if maintype == 'text':
                logger.debug('Attaching text file ' + fname)
                fp = open(fname)
                attach = MIMEText(fp.read(), _subtype=subtype)
                fp.close()
            elif maintype == 'image':
                logger.debug('Attaching image file ' + fname)
                fp = open(fname, 'rb')
                attach = MIMEImage(fp.read(), _subtype=subtype)
                fp.close()
            else:
                logger.debug('Attaching unknown file ' + fname)
                fp = open(fname, 'rb')
                attach = MIMEBase(maintype, subtype)
                attach.set_payload(fp.read())
                fp.close()
                # Encode the payload using Base64
                encoders.encode_base64(attach)
            # Set the filename parameter
            if a.get_desired_name() is None:
                filename = os.path.basename(fname)
            else:
                filename = a.get_desired_name()
            attach.add_header('Content-Disposition', 'attachment',
                              filename=filename)
            msg_root.attach(attach)

        return msg_root

    def _build_mime_message(self, from_address, to_list,
                            subject, message, reply_to, attachments):
        """Creates MIMEText object and returns it
        :param from_address: from email address
        :param to_list: list of email addresses as strings to send email to
        :param subject: Subject of email
        :param message: Body of email
        :returns MIMEText object upon success or None upon error
        """
        msg_root = MIMEMultipart("alternative")

        msg_root['Subject'] = subject
        msg_root['From'] = from_address
        msg_root['To'] = ", ".join(to_list)
        if reply_to is not None:
            msg_root.add_header('reply-to', reply_to)

        msg_root.attach(MIMEText(message, "plain"))
        updated_msg_root = self._append_attachments(msg_root, attachments)
        updated_msg_root.attach(MIMEText('<pre>\n' + message + '\n</pre>\n',
                                         "html"))
        updated_msg_root.preamble = ('You need a MIME enabled mail reader '
                                     'to see this message')

        return updated_msg_root
