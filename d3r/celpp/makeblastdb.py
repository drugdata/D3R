# -*- coding: utf-8 -*-

__author__ = 'churas'

import os.path
import logging
import time
from d3r.celpp.task import D3RTask
from d3r.celpp import util

logger = logging.getLogger(__name__)


class MakeBlastDBTask(D3RTask):
    """Represents Blastable database

       Downloads pdb_seqres.txt from rcsb, gunzips it and
       runs NCBI makeblastdb to generate a blast database
       """

    PDB_SEQRES_TXT = "pdb_seqres.txt"
    PDB_SEQRES_TXT_GZ = (PDB_SEQRES_TXT + ".gz")

    # PDBID codes are currently 4 characters so plus 1 is 5
    LENGTH_OF_PDBID_PLUS_ONE = 5

    TASK_NAME = 'makeblastdb'

    def __init__(self, path, args):
        """Constructor

           Constructor sets name to makeblastdb and stage is set to 1
        """
        super(MakeBlastDBTask, self).__init__(path, args)
        self.set_name(MakeBlastDBTask.TASK_NAME)
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._maxretries = 3
        self._retrysleep = 1

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files in addition to stderr/stdout files
           if they are found on the filesystem

           pdb_seqres.txt.gz

           :returns: list of files that can be uploaded.
        """
        # get the stderr/stdout files
        file_list = super(MakeBlastDBTask, self).get_uploadable_files()

        pdbfile = self.get_pdb_seqres_txt_gz()

        if os.path.isfile(pdbfile):
            file_list.append(pdbfile)

        return file_list

    def get_pdb_seqres_txt_gz(self):
        """Returns path to pdb_seqres.txt.gz file
        :return: full path to MakeBlastDBTask.PDB_SEQRES_TXT_GZ file
        """
        return os.path.join(self.get_dir(), MakeBlastDBTask.PDB_SEQRES_TXT_GZ)

    def get_pdb_seqres_txt(self):
        """Returns path to pdb_seqres.txt file
        :return: full path to MakeBlastDBTask.PDB_SEQRES_TXT file
        """
        return os.path.join(self.get_dir(), MakeBlastDBTask.PDB_SEQRES_TXT)

    def _get_sequence_count_message(self):
        """Returns number of fasta sequences in get_pdb_seqres_txt() file

           Examines file set by get_pdb_seqres_txt and counts number of
           lines starting with > character which denotes the id of a fasta
           sequence entry
           :return: string of format: # sequence(s): XXX where is
                    a number upon success or the message Error unable to parse
                    file upon error
        """
        logger.debug('Examing ' + self.get_pdb_seqres_txt() +
                     ' to get sequence count')
        seq_count = 0
        try:
            f = open(self.get_pdb_seqres_txt(), 'r')
            for line in f:
                if line[0] == '>':
                    seq_count += 1
            f.close()
            return '# sequence(s): ' + str(seq_count)
        except:
            logger.exception('Caught exception trying to get sequence count')
        return '# sequence(s): Error unable to parse file'

    def get_set_of_pbdid_from_pdb_seqres_txt(self):
        """Parses MakeBlastDBTask.PDB_SEQRES_TXT to get set of PBDIDs

           This method parses `MakeBlastDBTask.PDB_SEQRES_TXT` fasta file to
           get a unique set of PBDIDs from the deflines (lines starting with >)

           Format of defline:

           >PBDID_# ###:##### length:###  ########

           The # above changes and PBDID is the PBDID code.

           This method will grab the PBDID which is 4 characters long,
           but for safety will parse on _
           character and return those PBDID as a set.  The PDBIDs will
           be made uppercase with upper() string function

           :returns: set of PBDIDS from sequence file where PDBIDs set
           to upper case via upper()
        """
        pdbid_set = set()
        start_time = int(time.time())
        logger.debug('Examing ' + self.get_pdb_seqres_txt() +
                     ' to get set of PDBIDs')

        if not os.path.isfile(self.get_pdb_seqres_txt()):
            logger.warning('')
            return set()
        try:
            f = open(self.get_pdb_seqres_txt(), 'r')
            for line in f:
                if line[0] == '>':
                    underscore_pos = line.find('_')
                    if underscore_pos != \
                            MakeBlastDBTask.LENGTH_OF_PDBID_PLUS_ONE:
                        logger.error('Skipping entry... defline in sequence '
                                     'file missing '
                                     'underscore at position 5: ' + line)
                        continue

                    pdbid_set.add(line[1:underscore_pos].upper())
            f.close()
            logger.debug('Found ' + str(len(pdbid_set)) +
                         ' PDBIDs from sequence.  Parsing took ' +
                         str(int(time.time()) - start_time) +
                         ' seconds.')
            return pdbid_set
        except:
            logger.exception('Caught exception trying to get set of PDBIDs')
            return set()
        return pdbid_set

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
        super(MakeBlastDBTask, self).run()

        if self._can_run is False:
            logger.debug(self.get_dir_name() +
                         ' cannot run cause _can_run flag ' +
                         'is False')
            return

        try:
            logger.debug('pdbsequrl set to ' + self.get_args().pdbsequrl)
        except AttributeError:
            self.set_error('cannot download files cause pdbsequrl not set')
            self.end()
            return

        try:
            logger.debug('makeblastdb set to ' + self.get_args().makeblastdb)
        except AttributeError:
            self.set_error('cannot make blast database cause '
                           'makeblastdb not set')
            self.end()
            return

        download_path = self.get_pdb_seqres_txt_gz()
        url = self.get_args().pdbsequrl
        try:
            util.download_url_to_file(url,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

        except Exception:
            logger.exception('Caught Exception trying to download file: ' +
                             url)
            self.set_error('Unable to download file: ' + url)
            self.end()
            return

        try:
            # call util method to decompress file
            util.gunzip_file(download_path, self.get_pdb_seqres_txt())
        except IOError:
            logger.exception('Caught Exception trying to gunzip file: ' +
                             self.get_pdb_seqres_txt())
            self.set_error('Unable to uncompress file: ' +
                           self.get_pdb_seqres_txt())
            self.end()
            return

        # Run makeblastdb
        cmd_to_run = (self.get_args().makeblastdb + ' -in ' +
                      self.get_pdb_seqres_txt() +
                      ' -out ' + os.path.join(self.get_dir(), 'pdb_db') +
                      ' -dbtype prot')

        makeblastdb_name = os.path.basename(self.get_args().makeblastdb)

        self.run_external_command(makeblastdb_name, cmd_to_run, True)

        self.append_to_email_log('\n' + self._get_sequence_count_message() +
                                 '\n')
        # assess the result
        self.end()
