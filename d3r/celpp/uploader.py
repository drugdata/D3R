#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time
from ftpretty import ftpretty

logger = logging.getLogger(__name__)


class FileUploader(object):
    """Defines interface for uploading files
    """

    def upload_files(self, list_of_files):
        """Uploads files in `list_of_files`

           This method SHOULD be overridden in subclass
           :raises: NotImplementedError
           :returns: True upon success, false otherwise.
        """
        logger.debug('Running dummy get_upload_summary method')
        return True

    def get_upload_summary(self):
        """Gets summary of previous `upload_files` invocation
            This method SHOULD be overridden in subclass
            :raises: NotImplementedError
            :returns: Human readable string summary
        """
        logger.debug('Running dummy get_upload_summary method')
        return 'Nothing to report since this does not do anything'


class FtpFileUploader(FileUploader):
    """Implements FileUploader class by enabling upload of files to ftp

    """

    def __init__(self, ftp_config):
        """Constructor

        """
        # TODO need to parse ftp_config file
        self._connect_timeout = 60
        self._ftp = None
        self._files_transferred = 0
        self._bytes_transferred = 0
        self._duration = 0
        if ftp_config is None:
            return
        # self._ftp_remote_dir = ftp_remote_dir
        # self._ftp_host = ftp_host
        # self._ftp_user = ftp_user
        # self._ftp_pass = ftp_pass


    def set_ftp_connection(self, ftp_con):
        """Lets caller use library other then ftpretty for ftp connections

           If set it is assumed `ftp_con` is already connected to
           ftp server and it supports the API defined by ftpretty:
           ftp_con.close()
           ftp_con.put('/pathtofile','/destinationpath')
        """
        self._alt_ftp_con = ftp_con

    def set_ftp_remote_dir(self, remote_dir):
        """Sets the remote directory where ftp files will be uploaded to
        """
        self._remote_dir = remote_dir

    def get_ftp_remote_dir(self):
        return self._remote_dir

    def set_ftp_host(self, ftp_host):
        """Sets ftp host to connect to
        :param: host should be hostname ie ftp.box.com
        """
        self._ftp_host = ftp_host

    def get_ftp_host(self):
        """Gets ftp host
        """
        return self._host

    def set_ftp_user(self, ftp_user):
        self._ftp_user = ftp_user

    def get_ftp_user(self):
        return self._ftp_user

    def set_ftp_password(self, ftp_pass):
        self._ftp_pass = ftp_pass

    def get_ftp_password(self):
        return self._ftp_pass

    def set_connect_timeout(self, timeout):
        self._connect_timeout = timeout

    def get_connect_timeout(self):
        return self._connect_timeout

    def _connect(self):
        if self._alt_ftp_con is None:
            self._ftp = ftpretty(self.get_ftp_host(),
                                 self.get_ftp_user(),
                                 self.get_ftp_password(),
                                 timeout=self.get_connect_timeout())
            return
        self._ftp = self._alt_ftp_con

    def _disconnect(self):
        if self._ftp is None:
            return

        try:
            self._ftp.close()
        except:
            logger.exception('Caught exception attempting to close connection')

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

        self._ftp.put(file, os.path.join(self.get_ftp_remote_dir(),
                                         file))

        self._bytes_transferred += os.path.getsize(file)
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

           For each file this method first verifies the
           parent directory on the remote server exists.
           if not, the method attempts to create the
           necessary directories before attempting an
           upload.

           :returns: True upon success, False otherwise
        """
        start_time = int(time.time())
        try:
            if list_of_files is None:
                logger.warning('list of files is None')
                return True

            if len(list_of_files) is 0:
                logger.debug('No files to upload')
                return True

            # Try to connect
            try:
                self._connect()
            except:
                logger.exception('Unable to connect to ftp host')
                return False

            try:
                for file in list_of_files:
                    self._upload_file(file)
            except:
                logger.exception('Caught exception')
                return False
            finally:
                self._disconnect()
            return True
        finally:
            self._duration = int(time.time())  - start_time

    def get_upload_summary(self):
        """Gets summary of previous `upload_files` invocation
            This method SHOULD be overridden in subclass
            :raises: NotImplementedError
            :returns: Human readable string summary
        """
        logger.debug('Running dummy get_upload_summary method')

        return (str(self._files_transferred) + ' (' +
                str(self._bytes_transferred) + ' bytes) files uploaded in ' +
                str(self._duration) + ' seconds to host ' +
                self.get_ftp_host() + ':' + self.get_ftp_remote_dir())
