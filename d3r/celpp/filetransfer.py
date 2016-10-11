#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time
from ftpretty import ftpretty
import easywebdav

logger = logging.getLogger(__name__)


class InvalidFtpConfigException(Exception):
    """Thrown when an invalid ftp configuration is passed to FtpFileUploader
    """
    pass


class FileTransfer(object):
    """FileTransfer
    """
    HOST = 'host'
    USER = 'user'
    PASS = 'pass'
    PATH = 'path'
    CONTESTANTID = 'contestantid'
    CHALLENGEPATH = 'challengepath'
    SUBMISSIONPATH = 'submissionpath'

    def __init__(self, config):
        """Constructor

        """
        self._connect_timeout = 60
        self._ftp = None
        self._files_transferred = 0
        self._bytes_transferred = 0
        self._duration = 0
        self._alt_ftp_con = None
        self._ftp_host = None
        self._ftp_user = None
        self._ftp_pass = None
        self._remote_dir = None
        self._error_msg = None
        self._contestant_id = None
        self._remote_challenge_dir = None
        self._remote_submission_dir = None

        if config is not None:
            self._parse_config(config)

    def _parse_config(self, config):
        """Parses ftp config file for user credentials

           The parsed values are stored internally in this
           class.

           Expected format of `ftp_config` file

           host <HOST>
           user <USERNAME>
           pass <PASSWORD>
           path <BASE PATH ON REMOTE SERVER ie /foo>
           contestantid <CELPP CONTESTANT ID>
           challengepath <BASE PATH ON REMOTE SERVER ie /challenge>
           submissionpath <BASE PATH ON REMOTE SERVER ie /submission>

           Example:

           host ftp.box.com
           user bob@bob.com
           pass 12345
           path /upload
           contestantid 24680
           challengepath /challenge
           submissionpath /usersubmissions

           The above format matches the standard used
           by NCFTP with the exception of `path` and `contestantid`
           which are custom to this class

           :param ftp_config: Path to ftp config file
           :raises IOError: If there was an error opening the file
        """
        if config is None:
            raise InvalidFtpConfigException('No ftp file specified')

        f = None
        try:
            f = open(config, 'r')
            for line in f:
                split_line = line.split(' ')
                if not len(split_line) == 2:
                    continue
                if split_line[0] == FileTransfer.HOST:
                    self.set_host(split_line[1].rstrip())
                elif split_line[0] == FileTransfer.USER:
                    self.set_user(split_line[1].rstrip())
                elif split_line[0] == FileTransfer.PASS:
                    self.set_password(split_line[1].rstrip())
                elif split_line[0] == FileTransfer.PATH:
                    self.set_remote_dir(split_line[1].rstrip())
                elif split_line[0] == FileTransfer.CONTESTANTID:
                    self.set_contestant_id(split_line[1].rstrip())
                elif split_line[0] == FileTransfer.CHALLENGEPATH:
                    self.set_remote_challenge_dir(split_line[1].rstrip())
                elif split_line[0] == FileTransfer.SUBMISSIONPATH:
                    self.set_remote_submission_dir(split_line[1].rstrip())

        finally:
            if f is not None:
                f.close()

    def set_connection(self, ftp_con):
        """Lets caller use library other then ftpretty for ftp connections

           If set it is assumed `ftp_con` is already connected to
           ftp server and it supports the API defined by ftpretty:
           ftp_con.close()
           ftp_con.put('/pathtofile','/destinationpath')
        """
        self._alt_ftp_con = ftp_con

    def get_contestant_id(self):
        """Gets the contestant id
        """
        return self._contestant_id

    def set_contestant_id(self, contestant_id):
        """Sets the contestant id
        :param contestant_id: Id of contestant any value passed in other
        then None will be converted to a string
        """
        if contestant_id is None:
            self._contestant_id = None
            return

        self._contestant_id = str(contestant_id)

    def set_remote_challenge_dir(self, challenge_dir):
        """Sets the remote challenge directory for ftp upload
        """
        self._remote_challenge_dir = challenge_dir

    def get_remote_challenge_dir(self):
        """Gets the remote challenge directory for ftp upload
        """
        return self._remote_challenge_dir

    def set_remote_submission_dir(self, submission_dir):
        """Sets the remote submission directory for ftp download
        """
        self._remote_submission_dir = submission_dir

    def get_remote_submission_dir(self):
        """Gets the remote submission directory for ftp download
        """
        return self._remote_submission_dir

    def set_remote_dir(self, remote_dir):
        """Sets the remote directory where ftp files will be uploaded to
        """
        self._remote_dir = remote_dir

    def get_remote_dir(self):
        """Gets remote directory prefix path
           :returns: Remote directory prefix path, if value is None empty str
                     is returned
        """
        if self._remote_dir is None:
            return ''
        return self._remote_dir

    def set_host(self, ftp_host):
        """Sets ftp host to connect to
        :param: host should be hostname ie ftp.box.com
        """
        self._ftp_host = ftp_host

    def get_host(self):
        """Gets ftp host
        """
        return self._ftp_host

    def set_user(self, ftp_user):
        self._ftp_user = ftp_user

    def get_user(self):
        return self._ftp_user

    def set_password(self, ftp_pass):
        self._ftp_pass = ftp_pass

    def get_password(self):
        return self._ftp_pass

    def set_connect_timeout(self, timeout):
        self._connect_timeout = timeout

    def get_connect_timeout(self):
        return self._connect_timeout

    def get_error_msg(self):
        return self._error_msg

    def connect(self):
        """Dummy implementation
        """
        return True

    def disconnect(self):
        """Dummy implementation
        """
        pass

    def delete_file(self, remote_file):
        """Dummy implementation
        """
        self._error_msg = 'delete_file not implemented'
        return False

    def download_file(self, remote_file, local_file):
        """Dummy implementation
        """
        self._error_msg = 'download_file not implemented'
        return False

    def list_dirs(self, remote_dir):
        """Dummy
        """
        self._error_msg = 'list_dirs not implemented'
        return None

    def list_files(self, remote_dir):
        """Dummy
        """
        self._error_msg = 'list_files not implemented'
        return None

    def upload_file_direct(self, file, remote_dir, remote_file_name):
        """Dummy
        """
        self._error_msg = 'upload_file_direct not implemented'
        return False

    def upload_files(self, list_of_files):
        """Dummy
        """
        self._error_msg = 'upload_files not implemented'
        return False

    def get_upload_summary(self):
        """Gets summary of previous `upload_files` invocation

            :returns: Human readable string summary of format
            # files (# bytes) files uploaded in # seconds to
            host HOST:REMOTE_DIR
        """
        summary = ''
        if self._error_msg is not None:
            summary = self._error_msg + '\n'

        try:
            logger.debug('host: ' + self.get_host())
            host = self.get_host()
        except TypeError:
            logger.exception('Host not set')
            host = 'Unset'

        remote_dir = self.get_remote_dir()
        logger.debug('remote dir: ' + remote_dir)

        summary += (str(self._files_transferred) + ' (' +
                    str(self._bytes_transferred) +
                    ' bytes) files uploaded in ' +
                    str(self._duration) + ' seconds to host ' +
                    host + ':' + remote_dir)
        logger.debug('upload summary: ' + summary)
        return summary


class FtpFileTransfer(FileTransfer):
    """Implements FileUploader class by enabling upload of files to ftp

    """
    def __init__(self, config):
        """Constructor
        """
        super(FtpFileTransfer, self).__init__(config)

    def connect(self):
        if self._alt_ftp_con is None:
            try:

                logger.debug('Connecting to ' +
                             str(self.get_host()) + ' with user ' +
                             str(self.get_user()))
                self._ftp = ftpretty(self.get_host(),
                                     self.get_user(),
                                     self.get_password(),
                                     timeout=self.get_connect_timeout())
                return True
            except:
                logger.exception('Unable to connect to ftp host')
                self._error_msg = 'Unable to connect to ftp host'
                return False
        self._ftp = self._alt_ftp_con
        return True

    def disconnect(self):
        if self._ftp is None:
            logger.debug('ftp connection is None, just returning')
            return

        if self._alt_ftp_con is not None:
            logger.debug('external ftp connection used, skipping close')
            return

        try:
            logger.debug('Attempting to close ftp connection')
            self._ftp.close()
        except:
            logger.exception('Caught exception attempting to close connection')

    def delete_file(self, remote_file):
        """Deletes file specified by `remote_file`

           This method will delete file specified by
           `remote_file`.  If there is an error information can be
           obtained by calling `self.get_error_msg()`  It is assumed
           `connect()` has been called on this object.
           :param remote_file: full path to remote file to delete
           :returns: True upon success, false otherwise
        """
        self._error_msg = None
        start_time = int(time.time())
        try:
            if remote_file is None:
                self._error_msg = 'remote_file None'
                return False

            remote_path = os.path.normpath(remote_file)
            logger.debug('Deleting file : ' + remote_path)

            try:
                result = self._ftp.delete(remote_path)
                logger.debug(' Delete operation returned: ' + str(result))
            except Exception as e:
                logger.exception('Caught exception deleting file' +
                                 remote_path)
                self._error_msg = 'Unable to delete ' + remote_path +\
                                  ' : ' + str(e)
                return False

            return True
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('Delete operation took ' +
                         str(self._duration) + ' seconds')

    def download_file(self, remote_file, local_file):
        """Downloads file specified by `remote_file` to `local_file`

           This method will download file specified by
           `remote_file` to `local_file.  If there is an error information can
           be obtained by calling `self.get_error_msg()`  It is assumed
           `connect()` has been called on this object.
           :param remote_file: full path to remote file to download
           :param local_file: full path to download remote file to
           :returns: True upon success, false otherwise
        """
        self._error_msg = None
        start_time = int(time.time())
        try:
            if remote_file is None:
                self._error_msg = 'remote_file None'
                return False

            if local_file is None:
                self._error_msg = 'local_file None'
                return False

            remote_path = os.path.normpath(remote_file)
            local_path = os.path.normpath(local_file)
            logger.debug('Downloading file : ' + remote_path + ' to ' +
                         local_path)

            try:
                self._ftp.get(remote_path, local=local_path)
            except Exception as e:
                logger.exception('Caught exception downloading file' +
                                 remote_path + ' to ' + local_path)
                self._error_msg = ('Unable to download ' + remote_path +
                                   ' to ' + local_path + ' : ' + str(e))
                return False

            return True
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('Download operation took ' +
                         str(self._duration) + ' seconds')

    def list_dirs(self, remote_dir):
        """Gets list of all directories found in `remote_dir`

           This method gets a list of directories in `remote_dir`
           If there is an error information can be
           obtained by calling `self.get_error_msg()`  It is assumed
           `connect()` has been called on this object.
           :param remote_dir: full path to remote directory to examine
           :returns: None for failure and list upon success
        """
        self._error_msg = None
        start_time = int(time.time())
        try:
            if remote_dir is None:
                self._error_msg = 'remote_dir None'
                return None

            remote_path = os.path.normpath(remote_dir)
            logger.debug('Examining : ' + remote_path)

            dirlist = []
            try:
                for paths in self._ftp.list(remote_dir, extra=True):
                    if paths['directory'] == 'd':
                        if paths['name'] == '.':
                            continue
                        if paths['name'] == '..':
                            continue
                        dirlist.append(paths['name'])

            except Exception as e:
                logger.exception('Caught exception examining ' +
                                 remote_path)
                self._error_msg = ('Unable to get directory list for ' +
                                   remote_path +
                                   ' : ' + str(e))
                return None

            return dirlist
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('list directory operation took ' +
                         str(self._duration) + ' seconds')

    def list_files(self, remote_dir):
        """Gets list of all files found in `remote_dir`

           This method gets a list of files in `remote_dir`
           If there is an error information can be
           obtained by calling `self.get_error_msg()`  It is assumed
           `connect()` has been called on this object.
           :param remote_dir: full path to remote directory to examine
           :returns: None for failure and list upon success
        """
        self._error_msg = None
        start_time = int(time.time())
        try:
            if remote_dir is None:
                self._error_msg = 'remote_dir None'
                return None

            remote_path = os.path.normpath(remote_dir)
            logger.debug('Examining : ' + remote_path)

            filelist = []
            try:
                for paths in self._ftp.list(remote_dir, extra=True):
                    if paths['directory'] == '-':
                        filelist.append(paths['name'])

            except Exception as e:
                logger.exception('Caught exception examining ' +
                                 remote_path)
                self._error_msg = ('Unable to get file list for ' +
                                   remote_path +
                                   ' : ' + str(e))
                return None

            return filelist
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('list file operation took ' +
                         str(self._duration) + ' seconds')

    def upload_file_direct(self, file, remote_dir, remote_file_name):
        """Uploads file to remote server

           This method will upload the file to the `remote_dir` using the
           `remote_file_name` as the file name.  If there is an error
           information can be obtained by calling `self.get_error_msg()`
           :param file: full path to file to upload
           :param remote_dir: full path to remote directory to upload file to
           :param remote_file_name: name to use for file uploaded
           :returns: True upon success, false otherwise
        """
        self._error_msg = None
        self._bytes_transferred = 0
        self._files_transferred = 0
        start_time = int(time.time())
        try:
            if file is None:
                self._error_msg = 'File passed in is None'
                return False

            if not os.path.isfile(file):
                self._error_msg = file + ' is not a file'
                return False

            if remote_dir is None:
                self._error_msg = 'remote_dir is None'
                return False

            if remote_file_name is None:
                self._error_msg = 'remote_file_name is None'
                return False

            remote_path = os.path.normpath(remote_dir +
                                           os.sep + remote_file_name)
            logger.debug('Direct uploading file: ' + file + ' to ' +
                         remote_path)

            try:
                size = self._ftp.put(file, remote_path)

                logger.debug('  Uploaded ' + str(size) + ' bytes')
                self._bytes_transferred += size
                self._files_transferred += 1
            except Exception as e:
                logger.exception('Caught exception direct uploading file' +
                                 file + ' to ' + remote_path)
                self._error_msg = 'Unable to upload ' + file + ' to ' +\
                                  remote_path + ' : ' + str(e)
                return False
            return True
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('End of direct_file_upload operation took ' +
                         str(self._duration) + ' seconds')

    def _upload_file(self, file):
        """Uploads file to remote server

           This method assumes a valid ftp connection
           has been made with self._ftp

           :param: file full path to file to upload
        """
        if file is None:
            logger.error('File passed in is None')
            return

        if not os.path.isfile(file):
            logger.warning(file + ' is not a file')
            return

        remote_path = os.path.normpath(self.get_remote_dir() +
                                       os.sep + file)
        logger.debug('Uploading file: ' + file + ' to ' + remote_path)

        size = self._ftp.put(file, remote_path)

        logger.debug('  Uploaded ' + str(size) + ' bytes')
        self._bytes_transferred += size
        self._files_transferred += 1

    def upload_files(self, list_of_files):
        """Uploads files in `list_of_files` to ftp server

           Using ftp configuration set in constructor this
           method connects to the ftp server and uploads
           the files in `list_of_files` in a mirror fashion.
           For example, if `remote_dir` is set to:

           /foo

           and the files in `list_of_files` is:

           /home/bob/data.1
           /home/joe/data.2
           /var/blah.txt

           The files would be uploaded as follows:

           /home/bob/data.1 => /foo/home/bob/data.1
           /home/joe/data.2 => /foo/home/joe/data.2
           /var/blah.txt => /foo/var/blah.txt

           ftpretty is used for ftp upload and it makes
           sure the directories exist on remote server
           before uploading.

           This method will return False as soon as an
           ftp upload raises an exception.

           :returns: True upon success, False otherwise
        """
        self._error_msg = None
        self._bytes_transferred = 0
        self._files_transferred = 0
        start_time = int(time.time())
        try:
            if list_of_files is None:
                logger.warning('list of files is None')
                self._error_msg = 'List of files passed in was None'
                return True

            if len(list_of_files) is 0:
                logger.debug('No files to upload')
                self._error_msg = 'No files to upload'
                return True

            try:
                logger.debug('Uploading ' + str(len(list_of_files)) + ' files')
                for file in list_of_files:
                    self._upload_file(file)
            except:
                logger.exception('Caught exception')
                self._error_msg = 'Error during upload'
                return False
            return True
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('End of upload_files operation took ' +
                         str(self._duration) + ' seconds')


class WebDavFileTransfer(FileTransfer):
    """Implements FileUploader class by enabling upload of files to WebDav

    """
    def __init__(self, config):
        """Constructor
        """
        super(WebDavFileTransfer, self).__init__(config)

    def connect(self):
        if self._alt_ftp_con is None:
            try:

                logger.debug('Connecting to ' +
                             str(self.get_host()) + ' with user ' +
                             str(self.get_user()))
                self._ftp = easywebdav.connect(self.get_host(),
                                               username=self.get_user(),
                                               password=self.get_password(),
                                               protocol='https')
                return True
            except:
                logger.exception('Unable to connect to host')
                self._error_msg = 'Unable to connect to host'
                return False
        self._ftp = self._alt_ftp_con
        return True

    def download_file(self, remote_file, local_file):
        """Downloads file specified by `remote_file` to `local_file`

           This method will download file specified by
           `remote_file` to `local_file.  If there is an error information can
           be obtained by calling `self.get_error_msg()`  It is assumed
           `connect()` has been called on this object.
           :param remote_file: full path to remote file to download
           :param local_file: full path to download remote file to
           :returns: True upon success, false otherwise
        """
        self._error_msg = None
        start_time = int(time.time())
        try:
            if remote_file is None:
                self._error_msg = 'remote_file None'
                return False

            if local_file is None:
                self._error_msg = 'local_file None'
                return False

            remote_path = os.path.normpath(remote_file)
            local_path = os.path.normpath(local_file)
            logger.debug('Downloading file : ' + remote_path + ' to ' +
                         local_path)

            try:
                self._ftp.download(remote_path, local_path)
            except Exception as e:
                logger.exception('Caught exception downloading file' +
                                 remote_path + ' to ' + local_path)
                self._error_msg = ('Unable to download ' + remote_path +
                                   ' to ' + local_path + ' : ' + str(e))
                return False

            return True
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('Download operation took ' +
                         str(self._duration) + ' seconds')

    def upload_file_direct(self, file, remote_dir, remote_file_name):
        """Uploads file to remote server

           This method will upload the file to the `remote_dir` using the
           `remote_file_name` as the file name.  If there is an error
           information can be obtained by calling `self.get_error_msg()`
           :param file: full path to file to upload
           :param remote_dir: full path to remote directory to upload file to
           :param remote_file_name: name to use for file uploaded
           :returns: True upon success, false otherwise
        """
        self._error_msg = None
        self._bytes_transferred = 0
        self._files_transferred = 0
        start_time = int(time.time())
        try:
            if file is None:
                self._error_msg = 'File passed in is None'
                return False

            if not os.path.isfile(file):
                self._error_msg = file + ' is not a file'
                return False

            if remote_dir is None:
                self._error_msg = 'remote_dir is None'
                return False

            if remote_file_name is None:
                self._error_msg = 'remote_file_name is None'
                return False

            remote_path = os.path.normpath(remote_dir +
                                           os.sep + remote_file_name)
            logger.debug('Direct uploading file: ' + file + ' to ' +
                         remote_path)

            try:
                self._ftp.upload(file, remote_path)
                size = os.path.getsize(file)
                logger.debug('  Uploaded ' + str(size) + ' bytes')
                self._bytes_transferred += size
                self._files_transferred += 1
            except Exception as e:
                logger.exception('Caught exception direct uploading file' +
                                 file + ' to ' + remote_path)
                self._error_msg = 'Unable to upload ' + file + ' to ' +\
                                  remote_path + ' : ' + str(e)
                return False
            return True
        finally:
            self._duration = int(time.time()) - start_time
            logger.debug('End of direct_file_upload operation took ' +
                         str(self._duration) + ' seconds')
