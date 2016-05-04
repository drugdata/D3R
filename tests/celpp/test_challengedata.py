__author__ = 'churas'

import unittest
import tempfile
import os.path
import re

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
                task._create_readme(os.path.join(temp_dir,'doesnotexist'))
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

                if re.match('^celpprunner version: Unknown.*$', line):
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

                if re.match('^celpprunner version: 1.0.*$', line):
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
            dimport = DataImportTask(temp_dir,params)
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

    """
    def test_tar_challenge_dir_no_final_log_and_no_candidates(self):
        self.assertEqual(1, 2)

    def test_tar_challenge_dir(self):
        self.assertEqual(1, 2)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
          self.assertEqual(1, 2)

        finally:
            shutil.rmtree(temp_dir)
    """
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

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.genchallenge = 'true'
            params.pdbdb = '/foo'
            params.version = '1'
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            chall = ChallengeDataTask(temp_dir, params)

            chall.run()
            self.assertEqual(chall.get_error(), None)
            # test files get created
            errfile = os.path.join(chall.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(chall.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(chall.get_dir(),
                                  'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(chall.get_dir(),
                                  'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
