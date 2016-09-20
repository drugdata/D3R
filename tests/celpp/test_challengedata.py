__author__ = 'churas'

import unittest
import tempfile
import stat
import re
import tarfile
from mock import Mock

from d3r.celpp.filetransfer import FtpFileTransfer

"""
test_proteinligprep
--------------------------------

Tests for `proteinligprep` module.
"""

import shutil
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.dataimport import DataImportTask


class TestChallengeDataTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test where directory doesn't even exist
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            self.assertEqual(task.get_uploadable_files(), [])

            # test on empty dir
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # test with tarfile
            tarfile = task.get_celpp_challenge_data_tar_file()

            open(tarfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(tarfile)

            # test with additional stderr/stdout files
            errfile = os.path.join(task.get_dir(),
                                   'genchallengedata.py.stderr')
            open(errfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(errfile)

            outfile = os.path.join(task.get_dir(),
                                   'genchallengedata.py.stdout')
            open(outfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(outfile)
            flist.index(errfile)
            flist.index(tarfile)

        finally:
            shutil.rmtree(temp_dir)

    def test_get_celpp_challenge_data_dir_name_path_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            self.assertEqual(task.get_celpp_challenge_data_dir_name(),
                             'celpp_week0_0')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_celpp_challenge_data_dir_name_not_standard_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            task.create_dir()
            self.assertEqual(task.get_celpp_challenge_data_dir_name(),
                             'celpp_week0_0')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_celpp_challenge_data_dir_name(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            cyear = os.path.join(temp_dir, '2016')
            cweek = os.path.join(cyear, 'dataset.week.5')
            os.mkdir(cyear)
            os.mkdir(cweek)
            task = ChallengeDataTask(cweek, params)
            task.create_dir()
            self.assertEqual(task.get_celpp_challenge_data_dir_name(),
                             'celpp_week5_2016')

        finally:
            shutil.rmtree(temp_dir)

    def test_create_challenge_dir_unable_to_make_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            cyear = os.path.join(temp_dir, '2016')
            cweek = os.path.join(cyear, 'dataset.week.5')
            os.mkdir(cyear)
            os.mkdir(cweek)
            task = ChallengeDataTask(cweek, params)
            task.create_dir()
            cdir = os.path.join(task.get_dir(),
                                task.get_celpp_challenge_data_dir_name())
            open(cdir, 'a').close()

            try:
                task._create_challenge_dir()
                self.fail('Expected OSError')
            except OSError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_create_challenge_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            cyear = os.path.join(temp_dir, '2016')
            cweek = os.path.join(cyear, 'dataset.week.5')
            os.mkdir(cyear)
            os.mkdir(cweek)
            task = ChallengeDataTask(cweek, params)
            task.create_dir()
            cdir = os.path.join(task.get_dir(),
                                task.get_celpp_challenge_data_dir_name())
            task._create_challenge_dir()
            self.assertEqual(os.path.isdir(cdir), True)
            task._create_challenge_dir()
            self.assertEqual(os.path.isdir(cdir), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_create_readme_unable_to_write_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = "1.0.0"
            task = ChallengeDataTask(temp_dir, params)
            try:
                task._create_readme(os.path.join(temp_dir, 'doesnotexist'))
                self.fail('Expected IOError')
            except IOError:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_create_readme_version_unset_and_no_summary_txt(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            task._create_readme(temp_dir)
            readme = os.path.join(temp_dir,
                                  ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)
            f = open(readme, 'r')
            found = False
            for line in f:

                if re.match('^Celpprunner version: Unknown.*$', line):
                    found = True
                    break

            f.close()
            self.assertEqual(found, True)

        finally:
            shutil.rmtree(temp_dir)

    def test_create_readme_version_empty_summary(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = '1.0'
            task = ChallengeDataTask(temp_dir, params)
            blast = BlastNFilterTask(temp_dir, params)
            blast.create_dir()
            sfile = blast.get_blastnfilter_summary_file()
            open(sfile, 'a').close()

            task._create_readme(temp_dir)
            readme = os.path.join(temp_dir,
                                  ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)
            f = open(readme, 'r')
            found = False
            for line in f:

                if re.match('^Celpprunner version: 1.0.*$', line):
                    found = True
                    break

            f.close()
            self.assertEqual(found, True)

        finally:
            shutil.rmtree(temp_dir)

    def test_create_readme(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = '1.0'
            task = ChallengeDataTask(temp_dir, params)
            blast = BlastNFilterTask(temp_dir, params)
            blast.create_dir()
            sfile = blast.get_blastnfilter_summary_file()
            f = open(sfile, 'w')
            f.write('hello there\n')
            f.flush()
            f.close()

            task._create_readme(temp_dir)
            readme = os.path.join(temp_dir,
                                  ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)
            f = open(readme, 'r')
            found = False
            for line in f:
                print line
                if re.match('^hello there.*$', line):
                    found = True
                    break

            f.close()
            self.assertEqual(found, True)

        finally:
            shutil.rmtree(temp_dir)

    def test_copy_over_tsv_files_no_dataimport(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            task.create_dir()
            challenge_dir = task._create_challenge_dir()
            self.assertEqual(os.path.isdir(challenge_dir), True)
            task._copy_over_tsv_files(challenge_dir)

        finally:
            shutil.rmtree(temp_dir)

    def test_copy_over_tsv_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            dimport = DataImportTask(temp_dir, params)
            dimport.create_dir()

            ctsv = dimport.get_crystalph_tsv()
            f = open(ctsv, 'w')
            f.write('crystal')
            f.flush()
            f.close()

            nonpoly = dimport.get_nonpolymer_tsv()
            f = open(nonpoly, 'w')
            f.write('nonpoly')
            f.flush()
            f.close()

            seq = dimport.get_sequence_tsv()
            f = open(seq, 'w')
            f.write('seq')
            f.flush()
            f.close()

            task = ChallengeDataTask(temp_dir, params)
            task.create_dir()
            challenge_dir = task._create_challenge_dir()
            self.assertEqual(os.path.isdir(challenge_dir), True)
            task._copy_over_tsv_files(challenge_dir)

            cop_ctsv = os.path.join(challenge_dir,
                                    DataImportTask.CRYSTALPH_TSV)
            self.assertEqual(os.path.isfile(cop_ctsv), True)
            f = open(cop_ctsv)
            self.assertEqual(f.readline(), 'crystal')
            f.close()

            cop_nonpoly = os.path.join(challenge_dir,
                                       DataImportTask.NONPOLYMER_TSV)
            self.assertEqual(os.path.isfile(cop_nonpoly), True)
            f = open(cop_nonpoly)
            self.assertEqual(f.readline(), 'nonpoly')
            f.close()

            cop_seq = os.path.join(challenge_dir,
                                   DataImportTask.SEQUENCE_TSV)
            self.assertEqual(os.path.isfile(cop_seq), True)
            f = open(cop_seq)
            self.assertEqual(f.readline(), 'seq')
            f.close()

        finally:
            shutil.rmtree(temp_dir)

    def test_tar_challenge_dir_unable_to_write_tarfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            task = ChallengeDataTask(temp_dir, params)

            chall_dir = task.get_celpp_challenge_data_dir_name()
            try:
                task._tar_challenge_dir(chall_dir)
                self.fail('Expected IOError')
            except IOError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_tar_challenge_dir_nothing_but_readme(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            task.create_dir()
            chall_dir = task._create_challenge_dir()
            task._create_readme(chall_dir)

            self.assertEqual(os.path.isdir(chall_dir), True)
            name = task.get_celpp_challenge_data_dir_name()
            tfile = task._tar_challenge_dir(name)
            self.assertEqual(os.path.isfile(tfile), True)

            foodir = os.path.join(temp_dir, 'foo')
            tar = tarfile.open(tfile, 'r:*')
            tar.extractall(path=foodir)
            tar.close()
            cdir = os.path.join(foodir, name)
            readme = os.path.join(cdir, ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)

        finally:
            shutil.rmtree(temp_dir)

    def make_fake_candidate_dir(self, basedir, targetid,
                                candidateid, candidateligandid):
        """creates fake candidate dir
        """
        namedir = os.path.join(basedir, targetid)
        os.mkdir(namedir)
        file_list = []
        fasta = os.path.join(namedir, targetid + '.fasta')
        self.write_fake_data_to_file(fasta)
        file_list.append(fasta.replace(basedir + '/', ''))
        txt = os.path.join(namedir, targetid + '.txt')
        self.write_fake_data_to_file(txt)
        file_list.append(txt.replace(basedir + '/', ''))
        largest = os.path.join(namedir,
                               self.generate_pdb_file_name('largest',
                                                           targetid,
                                                           candidateid,
                                                           candidateligandid))
        self.write_fake_data_to_file(largest)
        file_list.append(largest.replace(basedir + '/', ''))
        smallest = os.path.join(namedir,
                                self.generate_pdb_file_name('smallest',
                                                            targetid,
                                                            candidateid,
                                                            candidateligandid))
        self.write_fake_data_to_file(smallest)
        file_list.append(smallest.replace(basedir + '/', ''))
        apo = os.path.join(namedir,
                           self.generate_pdb_file_name('apo',
                                                       targetid,
                                                       candidateid,
                                                       candidateligandid))
        self.write_fake_data_to_file(apo)
        file_list.append(apo.replace(basedir + '/', ''))

        holo = os.path.join(namedir,
                            self.generate_pdb_file_name('holo',
                                                        targetid,
                                                        candidateid,
                                                        candidateligandid))
        self.write_fake_data_to_file(holo)
        file_list.append(holo.replace(basedir + '/', ''))

        lig_smi = os.path.join(namedir, 'lig_' + candidateligandid +
                               '.smi')
        self.write_fake_data_to_file(lig_smi)
        file_list.append(lig_smi.replace(basedir + '/', ''))

        lig_inchi = os.path.join(namedir, 'lig_' + candidateligandid +
                                 '.inchi')
        self.write_fake_data_to_file(lig_inchi)
        file_list.append(lig_inchi.replace(basedir + '/', ''))

        lig_mol = os.path.join(namedir, 'lig_' + candidateligandid +
                               '.mol')
        self.write_fake_data_to_file(lig_mol)
        file_list.append(lig_mol.replace(basedir + '/', ''))
        return file_list

    def write_fake_data_to_file(self, filepath):
        """Writes fake data into file passed in

        The data content is the name of the file
        as extracted from os.path.basename(filepath)
        :param filepath: path to file
        """
        f = open(filepath, 'w')
        f.write(os.path.basename(filepath))
        f.flush()
        f.close()

    def generate_pdb_file_name(self, type, targetid, candidateid,
                               candidateligandid):
        """Generates a pdb file name matching challenge data dir format
        """
        return (type + '-' + targetid + '_' + candidateid + '_' +
                candidateligandid + '.pdb')

    def test_tar_challenge_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            dimport = DataImportTask(temp_dir, params)
            dimport.create_dir()
            open(dimport.get_sequence_tsv(), 'a').close()
            open(dimport.get_nonpolymer_tsv(), 'a').close()
            open(dimport.get_crystalph_tsv(), 'a').close()

            task = ChallengeDataTask(temp_dir, params)
            task.create_dir()
            chall_dir = task._create_challenge_dir()

            final_log = os.path.join(chall_dir, 'final.log')
            open(final_log, 'a').close()

            task._create_readme(chall_dir)
            task._copy_over_tsv_files(chall_dir)

            # make a fake candidate
            file_list = self.make_fake_candidate_dir(chall_dir, '5hib',
                                                     '2eb2', 'CSX')

            self.assertEqual(os.path.isdir(chall_dir), True)
            name = task.get_celpp_challenge_data_dir_name()
            tfile = task._tar_challenge_dir(name)
            self.assertEqual(os.path.isfile(tfile), True)

            foodir = os.path.join(temp_dir, 'foo')
            tar = tarfile.open(tfile, 'r:*')
            tar.extractall(path=foodir)
            tar.close()
            cdir = os.path.join(foodir, name)
            readme = os.path.join(cdir, ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)

            final_log = os.path.join(foodir, 'final.log')
            self.assertEqual(os.path.isfile(final_log), False)

            for fname in file_list:
                chk = os.path.join(cdir, fname)
                self.assertEqual(os.path.isfile(chk), True, chk)

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no blast task found so it cannot run
            params = D3RParameters()
            task = ChallengeDataTask(temp_dir, params)
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             'blastnfilter task has notfound status')

            # blastn filter running
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.START_FILE),
                 'a').close()
            task = ChallengeDataTask(temp_dir, params)
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             'blastnfilter task has start status')

            # blastnfilter failed
            error_file = os.path.join(blastnfilter.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            task = ChallengeDataTask(temp_dir, params)
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             'blastnfilter task has error status')

            # blastnfilter success
            os.remove(error_file)
            open(os.path.join(blastnfilter.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            task = ChallengeDataTask(temp_dir, params)
            self.assertEqual(task.can_run(), True)
            self.assertEqual(task.get_error(), None)

            # proteinligprep task exists already
            task = ChallengeDataTask(temp_dir, params)
            task.create_dir()
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             task.get_dir_name() +
                             ' already exists and status is unknown')

            # proteinlibprep already complete
            task = ChallengeDataTask(temp_dir, params)
            open(os.path.join(task.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = '1'
            # return immediately cause can_run is false
            chall = ChallengeDataTask(temp_dir, params)
            chall.run()
            self.assertEqual(chall.get_error(),
                             'blastnfilter task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_genchallenge_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = '1'
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            chall = ChallengeDataTask(temp_dir, params)
            chall.run()
            self.assertEqual(chall.get_error(),
                             'genchallenge not set')
            # test files get created
            self.assertEqual(os.path.isdir(chall.get_dir()),
                             True)
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_pdbdb_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.genchallenge = 'false'
            params.version = '1'
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            chall = ChallengeDataTask(temp_dir, params)
            chall.run()
            self.assertEqual(chall.get_error(),
                             'pdbdb not set')
            # test files get created
            self.assertEqual(os.path.isdir(chall.get_dir()),
                             True)
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_genchallenge_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.genchallenge = 'false'
            params.pdbdb = '/foo'
            params.version = '1'
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            chall = ChallengeDataTask(temp_dir, params)

            chall.run()
            self.assertEqual(chall.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(chall.get_dir(),
                                  'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(chall.get_dir(),
                                  'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_genchallenge_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.genchallenge = '/bin/doesnotexist'
            params.pdbdb = '/foo'
            params.version = '1'
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            chall = ChallengeDataTask(temp_dir, params)

            chall.run()
            self.assertEqual(chall.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --candidatedir ' +
                             blastnfilter.get_dir() + ' --pdbdb ' +
                             '/foo --outdir ' +
                             chall.get_dir() +
                             '/' + chall.get_celpp_challenge_data_dir_name() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_challenge_file_challenge_file_is_none(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()

            chall = ChallengeDataTask(temp_dir, params)
            chall._upload_challenge_file(None)
            chall.set_file_transfer(params)
            self.assertEqual(chall.get_email_log(), 'challenge_file is None in'
                                                    ' _upload_challenge_file'
                                                    '\n')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_challenge_file_uploader_is_none(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            chall = ChallengeDataTask(temp_dir, params)
            chall._upload_challenge_file('/foo')
            self.assertEqual(chall.get_email_log(), 'No uploader available '
                                                    'to upload challenge '
                                                    'data\n')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_challenge_file_no_remote_challenge_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            chall = ChallengeDataTask(temp_dir, params)
            ftp = FtpFileTransfer(None)
            chall.set_file_transfer(ftp)
            chall._upload_challenge_file('/foo')
            self.assertEqual(chall.get_email_log(), 'No remote challenge'
                                                    ' directory set for '
                                                    'ftp upload\n')
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_challenge_file_uploader_upload_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            yeardir = os.path.join(temp_dir, '2016')
            os.mkdir(yeardir)
            weekdir = os.path.join(yeardir, 'dataset.week.50')
            os.mkdir(weekdir)
            chall = ChallengeDataTask(weekdir, params)
            chall.create_dir()

            mockftp = D3RParameters()
            mockftp.put = Mock(side_effect=IOError('hi'))
            ftp = FtpFileTransfer(None)
            ftp.set_remote_challenge_dir('/challenge')
            ftp.set_connection(mockftp)
            ftp.connect()
            chall.set_file_transfer(ftp)
            tarball = os.path.join(chall.get_dir(), 'celppweek50_2016.tar.gz')
            f = open(tarball, 'w')
            f.write('hi')
            f.flush()
            f.close()
            try:
                chall._upload_challenge_file(tarball)
                self.fail('Expected exception')
            except Exception as e:
                self.assertEqual(str(e),
                                 'Unable to upload ' + tarball +
                                 ' to /challenge/celppweek50_2016.tar.gz : ' +
                                 'hi')
            ftp.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_challenge_file_uploader_tar_succeeds_latest_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            yeardir = os.path.join(temp_dir, '2016')
            os.mkdir(yeardir)
            weekdir = os.path.join(yeardir, 'dataset.week.50')
            os.mkdir(weekdir)
            chall = ChallengeDataTask(weekdir, params)
            chall.create_dir()

            mockftp = D3RParameters()
            mockftp.put = Mock(side_effect=[3, IOError('hi')])
            ftp = FtpFileTransfer(None)
            ftp.set_remote_challenge_dir('/challenge')
            ftp.set_connection(mockftp)
            ftp.connect()
            chall.set_file_transfer(ftp)
            tarball = os.path.join(chall.get_dir(), 'celppweek50_2016.tar.gz')
            f = open(tarball, 'w')
            f.write('hi')
            f.flush()
            f.close()
            latest_file = os.path.join(chall.get_dir(),
                                       ChallengeDataTask.LATEST_TXT)
            try:
                chall._upload_challenge_file(tarball)
                self.fail('Expected exception')
            except Exception as e:
                self.assertEqual(str(e),
                                 'Unable to upload ' + latest_file +
                                 ' to /challenge/latest.txt : ' +
                                 'hi')

            f = open(latest_file, 'r')
            line = f.readline()
            self.assertEqual(line, 'celppweek50_2016.tar.gz')
            f.close()
            ftp.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def test_upload_challenge_file_uploader_successful(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            yeardir = os.path.join(temp_dir, '2017')
            os.mkdir(yeardir)
            weekdir = os.path.join(yeardir, 'dataset.week.1')
            os.mkdir(weekdir)
            chall = ChallengeDataTask(weekdir, params)
            chall.create_dir()

            mockftp = D3RParameters()
            mockftp.put = Mock(side_effect=[3, 5])
            ftp = FtpFileTransfer(None)
            ftp.set_remote_challenge_dir('/challenge')
            ftp.set_connection(mockftp)
            ftp.connect()
            chall.set_file_transfer(ftp)
            tarball = os.path.join(chall.get_dir(), 'celppweek1_2017.tar.gz')
            f = open(tarball, 'w')
            f.write('hi')
            f.flush()
            f.close()
            latest_file = os.path.join(chall.get_dir(),
                                       ChallengeDataTask.LATEST_TXT)
            chall._upload_challenge_file(tarball)

            f = open(latest_file, 'r')
            line = f.readline()
            self.assertEqual(line, 'celppweek1_2017.tar.gz')
            f.close()
            ftp.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def create_gen_challenge_script(self, temp_dir):
        script = os.path.join(temp_dir, 'genchallenge.py')
        f = open(script, 'w')
        f.write("""#! /usr/bin/env python
import argparse
import sys
import os

class D3RParameters(object):
    pass

parsed_arguments = D3RParameters()
desc = 'hi'
help_formatter = argparse.RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=help_formatter)
parser.add_argument("--outdir", help='blah')
parser.add_argument("--candidatedir", help='blah')
parser.add_argument("--pdbdb", help='blah')

theargs = parser.parse_args(sys.argv[1:], namespace=parsed_arguments)

open(os.path.join(theargs.outdir,'final.log'), 'a').close()

thedir = os.path.join(theargs.outdir, 'error_container')
os.mkdir(thedir)
open(os.path.join(thedir,'foo.txt'),'a').close()

thedir = os.path.join(theargs.outdir, '5hib')
os.mkdir(thedir)

open(os.path.join(thedir,'5hib.fasta'),'a').close()
open(os.path.join(thedir,'5hib.txt'),'a').close()
open(os.path.join(thedir,'apo-5hib_2eb2.pdb'),'a').close()
open(os.path.join(thedir,'holo-5hib_5hg8-634.pdb'),'a').close()
open(os.path.join(thedir,'largest-5hib_5em6-5Q3.pdb'),'a').close()
open(os.path.join(thedir,'smallest-5hib_4wrg-CSX.pdb'),'a').close()
open(os.path.join(thedir,'lig_63M.inchi'),'a').close()
open(os.path.join(thedir,'lig_63M.mol'),'a').close()
open(os.path.join(thedir,'lig_63M.smi'),'a').close()

thedir = os.path.join(theargs.outdir, '5hic')
os.mkdir(thedir)
open(os.path.join(thedir,'5hic.fasta'),'a').close()
open(os.path.join(thedir,'5hic.txt'),'a').close()
open(os.path.join(thedir,'apo-5hic_2eb2.pdb'),'a').close()
open(os.path.join(thedir,'holo-5hic_5hg8-634.pdb'),'a').close()
open(os.path.join(thedir,'largest-5hic_5em6-5Q3.pdb'),'a').close()
open(os.path.join(thedir,'smallest-5hic_4wrg-CSX.pdb'),'a').close()
open(os.path.join(thedir,'lig_63N.inchi'),'a').close()
open(os.path.join(thedir,'lig_63N.mol'),'a').close()
open(os.path.join(thedir,'lig_63N.smi'),'a').close()

""")
        f.flush()
        f.close()
        os.chmod(script, stat.S_IRWXU)
        return script

    def test_run_fails_cause_ftp_upload_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            script = self.create_gen_challenge_script(temp_dir)
            params = D3RParameters()
            params.genchallenge = script
            params.pdbdb = '/foo'
            params.version = '1'

            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()

            chall = ChallengeDataTask(temp_dir, params)
            mockftp = D3RParameters()
            mockftp.put = Mock(side_effect=[3, IOError('hi')])
            ftp = FtpFileTransfer(None)
            ftp.set_remote_challenge_dir('/challenge')
            ftp.set_connection(mockftp)
            ftp.connect()
            chall.set_file_transfer(ftp)

            dimport = DataImportTask(temp_dir, params)
            dimport.create_dir()

            ctsv = dimport.get_crystalph_tsv()
            f = open(ctsv, 'w')
            f.write('crystal')
            f.flush()
            f.close()

            nonpoly = dimport.get_nonpolymer_tsv()
            f = open(nonpoly, 'w')
            f.write('nonpoly')
            f.flush()
            f.close()

            seq = dimport.get_sequence_tsv()
            f = open(seq, 'w')
            f.write('seq')
            f.flush()
            f.close()

            chall.run()
            latest_file = os.path.join(chall.get_dir(),
                                       ChallengeDataTask.LATEST_TXT)
            self.assertEqual(chall.get_error(), 'Caught exception ' +
                             'Unable to upload ' + latest_file + ' to ' +
                             '/challenge/' + ChallengeDataTask.LATEST_TXT +
                             ' : hi')
            # verify test files get created
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            compfile = os.path.join(chall.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), False)
            ftp.disconnect()
        finally:
            shutil.rmtree(temp_dir)

    def test_run_success_with_ftp_upload(self):
        temp_dir = tempfile.mkdtemp()
        try:
            script = self.create_gen_challenge_script(temp_dir)
            params = D3RParameters()
            params.genchallenge = script
            params.pdbdb = '/foo'
            params.version = '1'

            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()

            chall = ChallengeDataTask(temp_dir, params)
            mockftp = D3RParameters()
            mockftp.put = Mock(side_effect=[3, 5])
            ftp = FtpFileTransfer(None)
            ftp.set_remote_challenge_dir('/challenge')
            ftp.set_connection(mockftp)
            chall.set_file_transfer(ftp)

            dimport = DataImportTask(temp_dir, params)
            dimport.create_dir()

            ctsv = dimport.get_crystalph_tsv()
            f = open(ctsv, 'w')
            f.write('crystal')
            f.flush()
            f.close()

            nonpoly = dimport.get_nonpolymer_tsv()
            f = open(nonpoly, 'w')
            f.write('nonpoly')
            f.flush()
            f.close()

            seq = dimport.get_sequence_tsv()
            f = open(seq, 'w')
            f.write('seq')
            f.flush()
            f.close()

            chall.run()
            self.assertEqual(chall.get_error(), None)
            # verify test files get created
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(chall.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            ftp.disconnect()

        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds_no_ftp(self):
        temp_dir = tempfile.mkdtemp()
        try:
            script = self.create_gen_challenge_script(temp_dir)
            params = D3RParameters()
            params.genchallenge = script
            params.pdbdb = '/foo'
            params.version = '1'

            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            chall = ChallengeDataTask(temp_dir, params)

            dimport = DataImportTask(temp_dir, params)
            dimport.create_dir()

            ctsv = dimport.get_crystalph_tsv()
            f = open(ctsv, 'w')
            f.write('crystal')
            f.flush()
            f.close()

            nonpoly = dimport.get_nonpolymer_tsv()
            f = open(nonpoly, 'w')
            f.write('nonpoly')
            f.flush()
            f.close()

            seq = dimport.get_sequence_tsv()
            f = open(seq, 'w')
            f.write('seq')
            f.flush()
            f.close()

            chall.run()
            self.assertEqual(chall.get_error(), None)
            # verify test files get created
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(chall.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(chall.get_dir(),
                                  'genchallenge.py.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(chall.get_dir(),
                                  'genchallenge.py.stdout')
            self.assertEqual(os.path.isfile(stdout), True)

            # verify challenge directory is created and
            # filled with valid files
            chall_dir = os.path.join(chall.get_dir(),
                                     chall.get_celpp_challenge_data_dir_name())
            self.assertEqual(os.path.isdir(chall_dir), True)
            readme = os.path.join(chall_dir, ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)

            crystal = os.path.join(chall_dir, DataImportTask.CRYSTALPH_TSV)
            self.assertEqual(os.path.isfile(crystal), True)

            seq = os.path.join(chall_dir, DataImportTask.SEQUENCE_TSV)
            self.assertEqual(os.path.isfile(seq), True)

            nonpoly = os.path.join(chall_dir, DataImportTask.NONPOLYMER_TSV)
            self.assertEqual(os.path.isfile(nonpoly), True)

            fivehibdir = os.path.join(chall_dir, '5hib')
            fivehibtxt = os.path.join(fivehibdir, '5hib.txt')
            self.assertEqual(os.path.isfile(fivehibtxt), True)
            fivehibfas = os.path.join(fivehibdir, '5hib.fasta')
            self.assertEqual(os.path.isfile(fivehibfas), True)

            fivehicdir = os.path.join(chall_dir, '5hic')
            fivehictxt = os.path.join(fivehicdir, '5hic.txt')
            self.assertEqual(os.path.isfile(fivehictxt), True)
            fivehicfas = os.path.join(fivehicdir, '5hic.fasta')
            self.assertEqual(os.path.isfile(fivehicfas), True)

            # verify tarball is created
            tfile = chall.get_celpp_challenge_data_tar_file()
            self.assertEqual(os.path.isfile(tfile), True)

            foodir = os.path.join(temp_dir, 'foo')
            os.mkdir(foodir)

            tar = tarfile.open(tfile, 'r:*')
            tar.extractall(path=foodir)
            tar.close()
            name = chall.get_celpp_challenge_data_dir_name()
            cdir = os.path.join(foodir, name)
            readme = os.path.join(cdir, ChallengeDataTask.README_TXT_FILE)
            self.assertEqual(os.path.isfile(readme), True)

            final = os.path.join(foodir, name, 'final.log')
            self.assertEqual(os.path.isfile(final), False)

            e_con = os.path.join(foodir, name, 'error_container')
            self.assertEqual(os.path.isdir(e_con), False)

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
