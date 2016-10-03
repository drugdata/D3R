__author__ = 'churas'

import os
import logging
import tarfile
import shutil

from d3r.celpp import util
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.dataimport import DataImportTask

logger = logging.getLogger(__name__)


class ChallengeDataTask(D3RTask):
    """Generates ChallengeData files

    """

    TAR_EXCLUDE_DIRS = ['error_container']
    TAR_EXCLUDE_FILES = ['final.log']
    TAR_GZ_SUFFIX = ".tar.gz"
    README_TXT_FILE = "readme.txt"
    LATEST_TXT = "latest.txt"

    README_BODY = """CELPP Weekly Pose Prediction Challenge
======================================

Celpprunner version: {version}
Week: {week}
Year: {year}

This tar file contains the CELPP weekly pose prediction challenge dataset.

Within this readme.txt is a description of the data in this tar file as well as
a summary of the Blastnfilter run which generated these candidates.

Tsv files downloaded from
=========================

http://www.wwpdb.org/files/new_release_structure_sequence.tsv
http://www.wwpdb.org/files/new_release_structure_nonpolymer.tsv
http://www.wwpdb.org/files/new_release_crystallization_pH.tsv

Structure of data overview
==========================

This tar file contains a set of directories set to the name of Targets. Targets
are proteins which have primary sequence released, but not 3D coordinates.

Within each directory are a set of Candidates.  Candiates are proteins with
similar structure to the Target that also have known 3D coordinates which can
be used for pose prediction.

For more information visit:

https://github.com/drugdata/D3R
              or
https://drugdesigndata.org/about/celpp


Structure of data
=================

Below is a definition of the files and directories within this tar file:

[file or directory <text within denote values that change>]

  -- Definition


 [readme.txt]

     -- Description of data and output from celpp blastnfilter stage of
        processing.

 [new_release_crystallization_pH.tsv]
 [new_release_structure_nonpolymer.tsv]
 [new_release_structure_sequence.tsv]

     -- Tsv files downloaded from: http://www.wwpdb.org/files

 [<target id>]/
               [<target id>.txt]

                  -- Summary of Blastnfilter results for target protein
                     with PDBID.

               [LMCSS-<target id>_<candidate id>-<candidate ligand id>.pdb]

                  -- Candidate protein for docking which:
                      1) Passes the Blastnfilter criteria

                      2) Contains the Ligand with the largest maximum common
                         substructure (MCSS) to the Target Ligand.

                         Note:  If multiple proteins founded, the protein
                                with the highest resolution will be picked.

               [SMCSS-<target id>_<candidate id>-<candidate ligand id>.pdb]

                  -- Candidate protein for docking which:
                      1) Passes the Blastnfilter criteria.

                      2) Contains the Ligand with the smallest maximum common
                         substructure (MCSS) to the Target Ligand.

                         Note:  If multiple proteins founded, the protein
                                with the highest resolution will be picked.

               [hiResHolo-<target id>_<candidate id>-<candidate ligand id>.pdb]

                   -- Candidate protein for docking which:
                      1) Passes the Blastnfilter criteria.

                      2) Has the highest resolution among all holo proteins.

               [hiResApo-<target id>_<candidate id>-<candidate ligand id>.pdb]

                   -- Candidate protein for docking which:
                      1) Passes the Blastnfilter criteria.

                      2) Has the highest resolution among all apo proteins.

               [hiTanimoto-<target id>_<candidate id>-<candidate ligand id>\
                                                                          .pdb]

                   -- Candidate protein for docking which:
                      1) Passes the Blastnfilter criteria.

                      2) Contains the Ligand with the highest structural
                         similarity in Tanimoto score to the Target Ligand.

                         Note:  If multiple proteins founded, the protein
                                with the highest resolution will be picked.

               [LMCSS-<target id>_<candidate id>-<candidate ligand id> \
                                                                      -lig.pdb]

                    -- Contains the 3D coordinate of the atoms for the ligand
                       in the MaxMCSS candidate (LMCSS) protein.

               [lig_<candidate ligand id>.smi]

                   -- Canonical smile string of the Target Ligand which will be
                      used in later docking.

               [lig_<candidate ligand id>.inchi]

                   -- Inchi string of the Target Ligand.

               [lig_<candidate ligand id>.mol]

                   -- 2D structure of the Target Ligand.

Blastnfilter Summary
====================

"""

    def __init__(self, path, args):
        super(ChallengeDataTask, self).__init__(path, args)
        self.set_name('challengedata')

        # Make stage number one higher then BlastNFilter Stage
        blast = BlastNFilterTask(path, args)
        self.set_stage(blast.get_stage() + 1)

        self.set_status(D3RTask.UNKNOWN_STATUS)
        self._challenge_tarball_filename = None
        self._week_num = util.get_celpp_week_number_from_path(self.get_path())
        self._year = util.get_celpp_year_from_path(self.get_dir())

    def get_uploadable_files(self):
        """Returns list of files that can be uploaded to remote server

           List will contain these files in addition to stderr/stdout files
           if they are found on the filesystem


           :returns: list of files that can be uploaded.
        """
        # get the stderr/stdout files
        file_list = super(ChallengeDataTask, self).get_uploadable_files()
        tarfile = self.get_celpp_challenge_data_tar_file()
        if os.path.isfile(tarfile):
            file_list.append(tarfile)
        else:
            logger.warning('No tar file found!!!')

        return file_list

    def get_celpp_challenge_data_tar_file(self):
        """Returns path to challenge tar ball
        """
        return os.path.join(self.get_dir(),
                            self.get_celpp_challenge_data_dir_name() +
                            ".tar.gz")

    def get_celpp_challenge_data_dir_name(self):
        """Returns path to celpp challenge data directory name
        """
        return 'celpp_week' + str(self._week_num) + '_' + str(self._year)

    def _create_challenge_dir(self):
        """Creates the challenge directory
        """
        challenge_dir = os.path.join(self.get_dir(),
                                     self.get_celpp_challenge_data_dir_name())

        logger.debug('Creating directory: ' + challenge_dir)
        if not os.path.isdir(challenge_dir):
            os.mkdir(challenge_dir)
        return challenge_dir

    def _create_readme(self, path):
        """Creates readme.txt file for task

        """
        ver = 'Unknown'
        try:
            ver = self.get_args().version
        except AttributeError:
            logger.warning('Version unset using Unknown')

        f = open(os.path.join(path, ChallengeDataTask.README_TXT_FILE), 'w')
        f.write(ChallengeDataTask.README_BODY.format(
            version=ver,
            week=self._week_num,
            year=self._year))

        blast = BlastNFilterTask(self.get_path(), self.get_args())
        summary_file = blast.get_blastnfilter_summary_file()

        # append summary.txt file
        if os.path.isfile(summary_file):
            sumfile = open(summary_file, 'r')
            for line in sumfile:
                f.write(line)
            sumfile.close()

        f.flush()
        f.close()

    def _copy_over_tsv_files(self, challenge_dir):
        """Copies over tsv files from `DataImportTask`
        """
        dataimport = DataImportTask(self.get_path(), self.get_args())

        crystal_dest = os.path.join(challenge_dir,
                                    DataImportTask.CRYSTALPH_TSV)

        if os.path.isfile(dataimport.get_crystalph_tsv()):
            logger.debug('Copying over ' + dataimport.get_crystalph_tsv() +
                         'to ' + crystal_dest)
            shutil.copyfile(dataimport.get_crystalph_tsv(), crystal_dest)
        else:
            logger.warning(dataimport.get_crystalph_tsv() +
                           ' file does not exist')

        nonpoly_dest = os.path.join(challenge_dir,
                                    DataImportTask.NONPOLYMER_TSV)

        if os.path.isfile(dataimport.get_nonpolymer_tsv()):
            logger.debug('Copying over ' + dataimport.get_nonpolymer_tsv() +
                         'to ' + nonpoly_dest)
            shutil.copyfile(dataimport.get_nonpolymer_tsv(), nonpoly_dest)
        else:
            logger.warning(dataimport.get_nonpolymer_tsv() +
                           ' file does not exist')

        seq_dest = os.path.join(challenge_dir,
                                DataImportTask.SEQUENCE_TSV)

        if os.path.isfile(dataimport.get_sequence_tsv()):
            logger.debug('Copying over ' + dataimport.get_sequence_tsv() +
                         'to ' + seq_dest)
            shutil.copyfile(dataimport.get_sequence_tsv(), seq_dest)
        else:
            logger.warning(dataimport.get_sequence_tsv() +
                           ' file does not exist')

    def _tar_challenge_dir(self, challenge_dir_name):
        """Creates compressed tarball of challenge directory
        """

        tfile = os.path.join(self.get_dir(),
                             challenge_dir_name +
                             ChallengeDataTask.TAR_GZ_SUFFIX)

        logger.debug('Creating tar file: ' + tfile)

        tar = tarfile.open(tfile, 'w:gz')
        challenge_dir = os.path.join(self.get_dir(), challenge_dir_name)

        for entry in os.listdir(challenge_dir):
            fullpath = os.path.join(challenge_dir, entry)
            if os.path.isdir(fullpath):
                if entry in ChallengeDataTask.TAR_EXCLUDE_DIRS:
                    logger.debug('Skipping insertion of ' + entry +
                                 ' into tar file')
                    continue

                logger.debug('Adding directory to tarball: ' + entry)
                tar.add(fullpath, arcname=challenge_dir_name + '/' + entry)
                continue
            if os.path.isfile(fullpath):
                if entry in ChallengeDataTask.TAR_EXCLUDE_FILES:
                    logger.debug('Skipping insertion of ' + entry +
                                 ' into tar file')
                    continue
                logger.debug('Adding file to tarball: ' + entry)
                tar.add(fullpath, challenge_dir_name + '/' + entry)

        tar.close()
        tfile_size = str(os.path.getsize(tfile))
        logger.debug('Tarfile created and is: ' +
                     tfile_size + 'bytes in size')
        self.append_to_email_log('Tarfile: ' + tfile + ' (' +
                                 tfile_size + ' bytes) created.')
        return tfile

    def can_run(self):
        """Determines if task can actually run

           This method first verifies the `BlastNFilterTask` task
           has `D3RTask.COMPLETE_STATUS` for
           status.  The method then verifies a `ProteinLigPrepTask` does
           not already exist.  If above is not true then self.set_error()
           is set with information about the issue
           :return: True if can run otherwise False
        """
        self._can_run = False
        self._error = None
        # check blast
        blastnfilter = BlastNFilterTask(self._path, self._args)
        blastnfilter.update_status_from_filesystem()
        if blastnfilter.get_status() != D3RTask.COMPLETE_STATUS:
            logger.info('Cannot run ' + self.get_name() + 'task ' +
                        'because ' + blastnfilter.get_name() + 'task' +
                        'has a status of ' + blastnfilter.get_status())
            self.set_error(blastnfilter.get_name() + ' task has ' +
                           blastnfilter.get_status() + ' status')
            return False

        # check this task is not complete and does not exist

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

    def _upload_challenge_file(self, challenge_file):
        """Uploads challenge file to remote ftp server
        """
        if challenge_file is None:
            logger.error('challenge_file is None')
            self.append_to_email_log('challenge_file is None in '
                                     '_upload_challenge_file\n')
            return

        uploader = self.get_file_transfer()
        if uploader is None:
            logger.warning('No uploader available to upload challenge data')
            self.append_to_email_log('No uploader available to upload '
                                     'challenge data\n')
            return
        try:
            remote_dir = uploader.get_remote_challenge_dir()
            if remote_dir is None:
                logger.warning('No remote challenge directory set for '
                               'ftp upload')
                self.append_to_email_log('No remote challenge directory set ' +
                                         'for ftp upload\n')
                return
            logger.debug('Attempting to upload ' + challenge_file + ' to '
                         + remote_dir)

            logger.debug('Connecting to remote server to upload challenge data'
                         'package')
            # connect to remote server
            uploader.connect()

            if uploader.upload_file_direct(challenge_file, remote_dir,
                                           os.path.basename(challenge_file)) \
                    is False:
                raise Exception(uploader.get_error_msg())

            self.append_to_email_log('\nChallenge tarball uploaded: \n'
                                     + uploader.get_upload_summary() + '\n')

            # create latest.txt file and upload to ftp server
            self._upload_latest_file(uploader, challenge_file, remote_dir)
        finally:
            uploader.disconnect()

    def _upload_latest_file(self, uploader, challenge_file, remote_dir):
        """Creates and uploads latest.txt file containing name of tarfile
        """
        chall_name = os.path.basename(challenge_file)
        latest_file = os.path.join(self.get_dir(),
                                   ChallengeDataTask.LATEST_TXT)
        f = open(latest_file, 'w')
        f.write(chall_name)
        f.flush()
        f.close()
        if uploader.upload_file_direct(latest_file, remote_dir,
                                       ChallengeDataTask.LATEST_TXT) is False:
            raise Exception(uploader.get_error_msg())

    def run(self):
        """Runs ChallengeDataTask after verifying blastnfilter was good

           Method requires can_run() to be called before hand with
           successful outcome
           Otherwise method invokes D3RTask.start then this method
           creates a directory and invokes blastnfilter script.  Upon
           completion results are analyzed and success or error status
           is set appropriately and D3RTask.end is invoked
           """
        super(ChallengeDataTask, self).run()

        if self._can_run is False:
            logger.debug(
                self.get_dir_name() + ' cannot run cause _can_run flag '
                                      'is False')
            return

        try:
            logger.debug('genchallenge set to ' +
                         self.get_args().genchallenge)
        except AttributeError:
            self.set_error('genchallenge not set')
            self.end()
            return

        try:
            logger.debug('pdbdb set to ' +
                         self.get_args().pdbdb)
        except AttributeError:
            self.set_error('pdbdb not set')
            self.end()
            return

        # create challenge dir
        challenge_dir = ''
        try:
            challenge_dir = self._create_challenge_dir()
            logger.debug('Made challenge dir ' + challenge_dir)
        except OSError:
                logger.exception('Problem making challenge dir ' +
                                 challenge_dir)
                self.set_error('Unable to create ' + challenge_dir)
                self.end()

        blastnfilter = BlastNFilterTask(self.get_path(),
                                        self.get_args())
        #
        # genchallengedata.py --candidatedir <path to stage.3.blastnfilter> \
        # --pdbdb <path to pdb database> \
        # --outdir <path to stage.4.challengedata/celpp_week#_####/>
        cmd_to_run = (self.get_args().genchallenge + ' --candidatedir ' +
                      blastnfilter.get_dir() +
                      ' --pdbdb ' + self.get_args().pdbdb +
                      ' --outdir ' + challenge_dir)

        genchallenge_name = os.path.basename(self.get_args().genchallenge)
        self.run_external_command(genchallenge_name, cmd_to_run,
                                  True)

        try:
            # create readme.txt file
            self._create_readme(challenge_dir)

            # copy over tsv files
            self._copy_over_tsv_files(challenge_dir)

            # tar and compress up celpp_week#_#### directory omitting final.log
            tfile = self._tar_challenge_dir(os.path.basename(challenge_dir))

            # upload tar file to remote server
            self._upload_challenge_file(tfile)

        except Exception as e:
            logger.exception('Caught exception')
            self.set_error('Caught exception ' + str(e))
        # assess the result
        logger.debug('About to complete challengedata.py')
        self.end()
        logger.debug('Completed challengedata.py')
