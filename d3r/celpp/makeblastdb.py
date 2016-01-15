# -*- coding: utf-8 -*-

__author__ = 'churas'

import os
import logging
from d3r.celpp.task import D3RTask
from d3r.celpp import util

logger = logging.getLogger(__name__)


class MakeBlastDBTask(D3RTask):
    """Represents Blastable database

       This object does not follow the standard D3RTask
       structure.  This is because the blastable database
       is created externally in a legacy structure.

       """

    PDB_SEQRES_TXT = "pdb_seqres.txt"
    PDB_SEQRES_TXT_GZ = (PDB_SEQRES_TXT + ".gz")

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
        self._maxretries = 3
        self._retrysleep = 1

    def get_pdb_seqres_txt_gz(self):
        """Returns path to pdb_seqres.txt.gz file
        :return: full path to MakeBlastDBTask.PDB_SEQRES_TXT_GZ file
        """
        return os.path(self.get_dir(), MakeBlastDBTask.PDB_SEQRES_TXT_GZ)

    def get_pdb_seqres_txt(self):
        """Returns path to pdb_seqres.txt file
        :return: full path to MakeBlastDBTask.PDB_SEQRES_TXT file
        """
        return os.path(self.get_dir(), MakeBlastDBTask.PDB_SEQRES_TXT)

    def create_dir(self):
        """Creates path set in get_path()

           """
        the_path = os.path.join(self._path, self.get_dir_name())
        os.makedirs(the_path)
        return the_path

    def get_dir_name(self):
        """Will always return current
        :return: 'current' string
        """
        return 'current'

    def can_run(self):
        """Determines if task can run

           Verifies this task does not already exist.
           If task does exist and does not have
           `D3RTask.COMPLETE_STATUS then self.get_error()
           is set with information about the issue
           :return: True if task can be run False otherwise
        """
        self._can_run = False
        self._error = None

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
        """Generates pdb blast database

           First downloads pdb_seqres.txt.gz file from rcsb ftp,
           then gunzips file, and finally runs makeblastdb
           on file to generate blast database
        """
        self(MakeBlastDBTask, self).run()

        try:
            logger.debug('pdbsequrl set to ' + self._args.pdbsequrl)
        except AttributeError:
            self.set_error('cannot download files cause pdbsequrl not set')
            self.end()
            return

        if self._can_run is False:
            logger.debug(self.get_dir_name() +
                         ' cannot run cause _can_run flag ' +
                         'is False')
            return

        download_path = self.get_pdb_seqres_txt_gz()
        url = self._args.pdbsequrl
        try:
            util.download_url_to_file(url,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

        except Exception:
            logger.exception('Caught Exception trying to download file(s)')
            self.set_error('Unable to download file')
            self.end()
            return

        # call util method to decompress file
        util.gunzip_file(download_path, self.get_pdb_seqres_txt())

        # Run makeblastdb
        cmd_to_run = ('makeblastdb -in ' + MakeBlastDBTask.PDB_SEQRES_TXT +
                      ' -out pdb_db -dbtype prot')
        logger.debug('Need to run ' + cmd_to_run)
        self.end()
        return
