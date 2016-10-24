# -*- coding: utf-8 -*-

__author__ = 'churas'

import logging
import os
import time
from d3r.celpp.task import D3RTask
from d3r.celpp import util
from d3r.celpp.makeblastdb import MakeBlastDBTask

logger = logging.getLogger(__name__)


class ImportRetryCountExceededError(Exception):
    """Error raised when import retry count has been exceeded
    """
    pass


class DataImportTask(D3RTask):
    """Represents DataImport Task
    This task downloads 3 files named new_release_structure_nonpolymer.tsv,
    new_release_structure_sequence.tsv, and
    new_release_crystallization_pH.tsv from the web
    """
    PARTICIPANT_LIST_CSV = "participant_list.csv"
    NONPOLYMER_TSV = "new_release_structure_nonpolymer.tsv"
    OLDSEQUENCE_TSV = "new_release_structure_sequence.tsv"
    SEQUENCE_TSV = "new_release_structure_sequence_canonical.tsv"
    CRYSTALPH_TSV = "new_release_crystallization_pH.tsv"
    COMPINCHI_ICH = "Components-inchi.ich"

    # standard to append to NON_POLYMER_TSV file
    NONPOLYMER_TSV_STANDARD = "1FCZ    156     InChI=1S/C24H26O3/c1-23(2)13-" \
                              "14-24(3,4)20-15-18(10-11-19(20)23)21(25)12-7-" \
                              "16-5-8-17(9-6-16)22(26)27/h5-12,15H,13-14H2,1" \
                              "-4H3,(H,26,27)/b12-7+\n"

    # standard to append to CRYSTALPH_TSV
    CRYSTALPH_TSV_STANDARD = "1FCZ	7\n"

    # standard to append to SEQUENCE_TSV
    SEQUENCE_TSV_STANDARD = "1FCZ    1       SPQLEELITKVSKAHQETFPSLCQLGKYTTN" \
                            "SSADHRVQLDLGLWDKFSELATKCIIKIVEFAKRLPGFTGLSIADQI" \
                            "TLLKAACLDILMLRICTRYTPEQDTMTFSDGLTLNRTQMHNAGFGPL" \
                            "TDLVFAFAGQLLPLEMDDTETGLLSAICLICGDRMDLEEPEKVDKLQ" \
                            "EPLLEALRLYARRRRPSQPYMFPRMLMKITDLRGISTKGAERAITLK" \
                            "MEIPGPMPPLIREMLE\n"

    TASK_NAME = 'dataimport'

    def __init__(self, path, args):
        super(DataImportTask, self).__init__(path, args)
        self.set_name(DataImportTask.TASK_NAME)
        makeblast = MakeBlastDBTask('', args)
        self.set_stage(makeblast.get_stage() + 1)
        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._maxretries = 3
        self._retrysleep = 1

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files in addition to stderr/stdout files
           if they are found on the filesystem

           new_release_structure_nonpolymer.tsv
           new_release_structure_sequence.tsv
           new_release_crystallization_pH.tsv
           Components-inchi.ich

           :returns: list of files that can be uploaded.
        """
        # get the stderr/stdout files
        file_list = super(DataImportTask, self).get_uploadable_files()

        compinchi = self.get_components_inchi_file()

        if os.path.isfile(compinchi):
            file_list.append(compinchi)

        crystal = self.get_crystalph_tsv()

        if os.path.isfile(crystal):
            file_list.append(crystal)

        nonpoly = self.get_nonpolymer_tsv()

        if os.path.isfile(nonpoly):
            file_list.append(nonpoly)

        seq = self.get_sequence_tsv()

        if os.path.isfile(seq):
            file_list.append(seq)

        return file_list

    def get_nonpolymer_tsv(self):
        """Returns path to new_release_structure_nonpolymer.tsv file
        :return: full path to DataImportTask.NONPOLYMER_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.NONPOLYMER_TSV)

    def get_oldsequence_tsv(self):
        """Returns path to new_release_structure_sequence.tsv file
        :return: full path to DataImportTask.SEQUENCE_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.OLDSEQUENCE_TSV)

    def get_sequence_tsv(self):
        """Returns path to new_release_structure_sequence_canonical.tsv file
        :return: full path to DataImportTask.SEQUENCE_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.SEQUENCE_TSV)

    def get_participant_list_csv(self):
        """Returns path to participant_list.csv file
        :returns: full path to DataImportTask.PARTICIPANT_LIST_CSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.PARTICIPANT_LIST_CSV)

    def get_crystalph_tsv(self):
        """Returns path to new_release_crystallization_pH.tsv file
        :return: full path to DataImportTask.CRYSTALPH_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.CRYSTALPH_TSV)

    def get_components_inchi_file(self):
        return os.path.join(self.get_dir(),
                            DataImportTask.COMPINCHI_ICH)

    def get_set_of_pdbid_from_crystalph_tsv(self):
        """Gets set of PDBID by parsing `DataImportTask.CRYSTALPH_TSV`

           Parses `DataImportTask.CRYSTALPH_TSV` to build a set of PDBIDs
           from the first 4 characters on each line minus the first line.

           Format space between columns appears to be tab:

            PDB_ID	_exptl_crystal_grow.pH
            <PDBID>	<PH>


           Example file:
            PDB_ID	_exptl_crystal_grow.pH
            4RFR	7.5
            4X09	6.5
            4XET	6.2
            4XF1	6.2
            4XF3	6.2

           :returns: set of PDBID in uppercase
        """
        pdbid_set = set()
        start_time = int(time.time())
        logger.debug('Examing ' + self.get_crystalph_tsv() +
                     ' to get set of PDBIDs')

        if not os.path.isfile(self.get_crystalph_tsv()):
            logger.warning('')
            return set()
        try:

            f = open(self.get_crystalph_tsv(), 'r')

            headerline = f.readline()
            if not headerline.startswith('PDB_ID'):
                logger.warning('First line in ' + self.get_crystalph_tsv() +
                               'did NOT start with PDBID_ID ('
                               + headerline + ')')
            for line in f:
                    pdbid_split = line.split()
                    if len(pdbid_split) != 2:
                        logger.error('Skipping entry... whitespace'
                                     ' split resulted in:' +
                                     str(len(pdbid_split)) +
                                     ' instead of expected 2: ' +
                                     line)
                        continue

                    pdbid_set.add(pdbid_split[0].upper())
            f.close()
            logger.debug('Found ' + str(len(pdbid_set)) +
                         ' PDBIDs from crystalph.  Parsing took ' +
                         str(int(time.time()) - start_time) +
                         ' seconds.')
            return pdbid_set
        except:
            logger.exception('Caught exception trying to get set of PDBIDs')
        return set()

    def get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres(self):
        """Gets set of PDBIDs that are in both tsv and sequence file

           Examines `DataImportTask.CRYSTALPH_TSV` and
           `MakeBlastDBTask.PDB_SEQRES_TXT` and returns a set of PDBIDs
           that are in both files
           :returns: set of PDBIDs uppercase that are in both files above
        """
        make_blastdb = MakeBlastDBTask(self._path, self._args)

        if not os.path.isfile(make_blastdb.get_pdb_seqres_txt()):
            logger.warning('No ' + make_blastdb.get_pdb_seqres_txt() +
                           ' file found')
            return set()

        c_pdbid_set = self.get_set_of_pdbid_from_crystalph_tsv()

        if len(c_pdbid_set) == 0:
            logger.warning('No PDBIds found in ' + self.get_crystalph_tsv())
            return set()

        seq_pdbid_set = make_blastdb.get_set_of_pbdid_from_pdb_seqres_txt()

        if len(seq_pdbid_set) == 0:
            logger.warning('No PDBIds found in ' +
                           make_blastdb.get_pdb_seqres_txt())
            return set()

        common_pdbid = set()

        # iterate through tsv pdb ids and return any found in
        # sequence pdb id set
        for id in c_pdbid_set:
            if id in seq_pdbid_set:
                common_pdbid.add(id)

        logger.debug('Found ' + str(len(common_pdbid)) + ' PDBIDs in ' +
                     self.get_crystalph_tsv() + ' and ' +
                     make_blastdb.get_pdb_seqres_txt())

        return common_pdbid

    def append_standard_to_files(self):
        """Appends standard 1FCZ to tsv files as mapped below

           NONPOLYMER_TSV_STANDARD is appended to `get_nonpolymer_tsv()`
           CRYSTALPH_TSV_STANDARD is appended to `get_crystalph_tsv()`
           SEQUENCE_TSV_STANDARD is appended to `get_sequence_tsv()`
           SEQUENCE_TSV_STANDARD is appended to `get_oldsequence_tsv()`
        """
        logger.debug('Appending 1FCZ standard to tsv files')
        try:
            util.append_string_to_file(self.get_nonpolymer_tsv(),
                                       DataImportTask.NONPOLYMER_TSV_STANDARD)
            util.append_string_to_file(self.get_sequence_tsv(),
                                       DataImportTask.SEQUENCE_TSV_STANDARD)
            util.append_string_to_file(self.get_oldsequence_tsv(),
                                       DataImportTask.SEQUENCE_TSV_STANDARD)
            util.append_string_to_file(self.get_crystalph_tsv(),
                                       DataImportTask.CRYSTALPH_TSV_STANDARD)

            self.append_to_email_log('\nAppending 1FCZ standard to tsv'
                                     ' files\n')
            nonpoly_count = util.get_file_line_count(self.get_nonpolymer_tsv())
            seq_count = util.get_file_line_count(self.get_sequence_tsv())
            crystal_count = util.get_file_line_count(self.get_crystalph_tsv())
            inchi_count = util.get_file_line_count(self.
                                                   get_components_inchi_file())

            self.append_to_email_log('Line counts:\n')
            self.append_to_email_log(str(inchi_count) + ' ' +
                                     DataImportTask.COMPINCHI_ICH + '\n')

            self.append_to_email_log(str(nonpoly_count) + ' ' +
                                     DataImportTask.NONPOLYMER_TSV + '\n')

            self.append_to_email_log(str(seq_count) + ' ' +
                                     DataImportTask.SEQUENCE_TSV + '\n')
            self.append_to_email_log(str(seq_count) + ' ' +
                                     DataImportTask.OLDSEQUENCE_TSV + '\n')
            self.append_to_email_log(str(crystal_count) + ' ' +
                                     DataImportTask.CRYSTALPH_TSV + '\n')

        except IOError:
            logger.exception('Error appending standard to file')

    def _wait_for_url_to_be_updated(self, url):
        """Waits for `url` to have been updated since start of celpp week
        """
        try:
            if self.get_args().skipimportwait is True:
                logger.debug('Skiping wait for tsv files to be updated')
                return
        except Exception:
            logger.exception('Skiping wait for tsv files to be updated')
            return

        try:
            importsleep = int(self.get_args().importsleep)
        except Exception:
            logger.exception('importsleep was not set using 1 second')
            importsleep = 1

        val = util.has_url_been_updated_since_start_of_celpp_week(url)
        counter = 0
        while val is False:
            if counter > self.get_args().importretry:
                raise ImportRetryCountExceededError(
                    url + ' has not been updated after ' +
                    str(counter*importsleep) +
                    ' seconds')
            logger.debug('Try #' + str(counter) + ' ' + url +
                         ' has not been updated. Sleeping ' +
                         str(importsleep) + ' seconds')
            time.sleep(importsleep)

            val = util.has_url_been_updated_since_start_of_celpp_week(url)
            counter += 1

    def can_run(self):
        """Determines if task can actually run

           The method verifies a `DataImportTask` does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None

        # check this task is not complete and does not exist

        self.update_status_from_filesystem()
        if self.get_status() == D3RTask.COMPLETE_STATUS:
            logger.debug("No work needed for " + self.get_name() +
                         " task")
            return False

        make_blastdb = MakeBlastDBTask(self._path, self._args)
        make_blastdb.update_status_from_filesystem()
        if make_blastdb.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + ' task ' +
                        'because ' + make_blastdb.get_name() + ' task' +
                        'has a status of ' + make_blastdb.get_status())
            self.set_error(make_blastdb.get_name() + ' task has ' +
                           make_blastdb.get_status() + ' status')
            return False

        if self.get_status() != D3RTask.NOTFOUND_STATUS:
            logger.warning(self.get_name() + " task was already " +
                           "attempted, but there was a problem")
            self.set_error(self.get_dir_name() + ' already exists and ' +
                           'status is ' + self.get_status())
            return False
        self._can_run = True
        return True

    def _download_files(self, url):

        try:
            nonpolyurl = url + '/' + DataImportTask.NONPOLYMER_TSV

            self._wait_for_url_to_be_updated(nonpolyurl)

            download_path = self.get_nonpolymer_tsv()
            util.download_url_to_file(nonpolyurl,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            sequence = url + '/' + DataImportTask.SEQUENCE_TSV
            self._wait_for_url_to_be_updated(sequence)

            download_path = self.get_sequence_tsv()
            util.download_url_to_file(sequence,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            oldsequence = url + '/' + DataImportTask.OLDSEQUENCE_TSV
            self._wait_for_url_to_be_updated(oldsequence)

            download_path = self.get_oldsequence_tsv()
            util.download_url_to_file(oldsequence,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            crystal = url + '/' + DataImportTask.CRYSTALPH_TSV

            self._wait_for_url_to_be_updated(crystal)
            download_path = self.get_crystalph_tsv()
            util.download_url_to_file(crystal,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)

            url = self._args.compinchi
            download_path = self.get_components_inchi_file()
            util.download_url_to_file(url +
                                      '/' + DataImportTask.COMPINCHI_ICH,
                                      download_path,
                                      self._maxretries,
                                      self._retrysleep)
            return True
        except Exception:
            logger.exception('Caught Exception trying to download file(s)')

            self.set_error('Unable to download file from ' +
                           url + ' to ' + download_path)
            # assess the result
            self.end()
        return False

    def _download_participant_list_csv(self):
        """Downloads from ftp the DataImportTask.PARTICIPANT_LIST_CSV file
           If there is an error. A message is appended to email log, but its
           not considered a failure
        """
        localfile = self.get_participant_list_csv()

        try:
            ft = self.get_file_transfer()
            ft.connect()
            remotefile = os.path.join(ft.get_remote_dir(),
                                      DataImportTask.PARTICIPANT_LIST_CSV)
            logger.debug('Downloading ' + remotefile + ' to ' +
                         localfile)
            ft.download_file(remotefile,
                             localfile)
        except Exception as e:
            logger.warning('Caught exception, but not considering this fatal')
            self.append_to_email_log('\nWARNING: Unable to download ' +
                                     DataImportTask.PARTICIPANT_LIST_CSV +
                                     ' which means external users will NOT '
                                     'get evaluation email : ' + str(e) + '\n')
            return
        finally:
            try:
                ft.disconnect()
            except Exception as e:
                logger.warning('Caught exception trying to disconnect')
        if not os.path.isfile(localfile):
            self.append_to_email_log('\nWARNING: ' +
                                     DataImportTask.PARTICIPANT_LIST_CSV +
                                     ' not downloaded which means external '
                                     'users will NOT get evaluation email\n')

    def run(self):
        """Downloads Components-inchi.ich

           Downloads 3 files from url specified in `self._args.pdbfileurl`
           namely the ones mentioned in the constructor to get_dir() directory.
           Each download will be retried up to self._maxretries time
           which is set in constructor.  After which set_error()
           will be set with message if file is still unable to be downloaded.
           """
        super(DataImportTask, self).run()

        try:
            logger.debug('pdbfileurl set to ' + self._args.pdbfileurl)
        except AttributeError:
            self.set_error('cannot download files cause pdbfileurl not set')
            self.end()
            return

        try:
            logger.debug('compinchi set to ' + self._args.compinchi)
        except AttributeError:
            self.set_error('cannot download files cause compinchi not set')
            self.end()
            return

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        url = self._args.pdbfileurl

        # if False is returned then the download failed
        # just return since the method has already
        # called end() and set_error()
        if self._download_files(url) is False:
            return

        # download the participant_list.csv file
        # if its missing no biggy
        self._download_participant_list_csv()

        # Compare TSV files with pdb_seqres.txt file to see if there
        # are any duplicates issue #15
        try:
            pdbs = self.get_set_of_pdbid_in_crystalph_tsv_and_pdb_seqres()
            if len(pdbs) == 0:
                self.append_to_email_log('\nFound no entries in ' +
                                         DataImportTask.CRYSTALPH_TSV +
                                         'and ' +
                                         MakeBlastDBTask.PDB_SEQRES_TXT +
                                         ' files\n')
            else:
                self.append_to_email_log('\nWARNING: Found ' +
                                         str(len(pdbs)) +
                                         ' PDBIDs in both ' +
                                         DataImportTask.CRYSTALPH_TSV +
                                         ' and ' +
                                         MakeBlastDBTask.PDB_SEQRES_TXT +
                                         ' files\n')
        except Exception:
            logger.exception('Caught exception comparing tsv with '
                             'sequence')
            self.append_to_email_log('\nError unable to examine ' +
                                     DataImportTask.CRYSTALPH_TSV +
                                     ' and/or ' +
                                     MakeBlastDBTask.PDB_SEQRES_TXT +
                                     ' files\n')

        # Append internal standard
        self.append_standard_to_files()
        self.end()
