# -*- coding: utf-8 -*-

import os
import time
import logging
import smtplib
import platform
import mimetypes
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
        except:
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
        except:
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
        try:
            program_name = self.get_args().program
        except:
            logger.debug('args.program  was not set.  '
                         'using __file__')
        return program_name + ' ' + self._get_program_version()

    def _get_program_version(self):
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
        if self.get_file_transfer() is None:
            return
        try:
            uploadable_files = self.get_uploadable_files()
            self.get_file_transfer().connect()
            self.get_file_transfer().upload_files(uploadable_files)
            summary = self.get_file_transfer().get_upload_summary()
            self.append_to_email_log('\n' + summary + '\n')
        except:
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
        sfile.write(self._get_program_version())
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
            except:
                logger.exception('Caught exception')

        self._upload_task()

        if self.get_status() == D3RTask.ERROR_STATUS:
            try:
                err_file = os.path.join(self.get_dir(),
                                        D3RTask.ERROR_FILE)
                if not os.path.isfile(err_file):
                    open(err_file, 'a').close()
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
        try:
            logger.debug('Checking path ' + path_to_check)
        except:
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
            returncode, out, err = util.run_external_command(cmd_to_run)
        except Exception as e:
            logger.exception("Error caught exception")
            self.set_status(D3RTask.ERROR_STATUS)
            self.set_error("Caught Exception trying to run " +
                           cmd_to_run + " : " + str(e))
            self.end()
            return 1

        self.write_to_file(err, command_name + D3RTask.STDERR_SUFFIX)
        self.write_to_file(out, command_name + D3RTask.STDOUT_SUFFIX)

        if returncode != 0:
            if command_failure_is_fatal:
                self.set_status(D3RTask.ERROR_STATUS)
                self.set_error("Non zero exit code: " + str(returncode) +
                               " received. Standard out: " + out +
                               " Standard error: " + err)
            else:
                self.append_to_email_log("Although considered non fatal " +
                                         "for processing of stage a " +
                                         "non zero exit code: " +
                                         str(returncode) +
                                         "received. Standard out: " + out +
                                         " Standard error : " + err)
        return returncode


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


class SmtpEmailer(object):
    """Simple Wrapper class to send email via smtplib
    """

    def __init__(self, smtp_host='localhost', port=25):
        """Constructor
        :param smtp_host: host of stmp server
        :param port: port stmp server is on
        """
        self._smtp_host = smtp_host
        self._port = port
        self._alt_smtp_server = None

    def set_alternate_smtp_server(self, server):
        """Sets alternate smtp server to use
        """
        self._alt_smtp_server = server

    def send_email(self, from_address, to_list, subject, message,
                   reply_to=None, attachments=None):
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
            mime_msg = self._build_mime_message(from_address, to_list,
                                                subject, message, reply_to,
                                                attachments)
            server = self._get_server()
            server.sendmail(from_address, to_list, mime_msg.as_string())
        except smtplib.SMTPConnectError as e:
            logger.exception('Caught exception')
            raise EmailSendError('Unable to connect to smtp server ' + str(e))
        except Exception as e:
            logger.exception('Caught exception')
            raise EmailSendError('Caught exception ' + str(e))
        finally:
            if server is not None:
                server.quit()

    def generate_from_address_using_login_and_host(self):
        """Creates from email address from login and hostname
        of machine running this script.
        :returns: from email address as string
        """
        hostname = platform.node()
        if hostname is '':
            hostname = 'localhost'
        return os.getlogin() + '@' + hostname

    def _get_server(self):
        """Gets Smtp server
        :returns: smtplib.SMTP object unless set_alternate_smtp_server
        was set in
        which case that object is returned
        """
        if self._alt_smtp_server is not None:
            return self._alt_smtp_server

        return smtplib.SMTP(self._smtp_host,
                            self._port)

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
