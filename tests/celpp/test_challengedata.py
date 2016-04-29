__author__ = 'churas'

import unittest
import tempfile
import os.path

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


class TestChallengeDataTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        self.assertEqual(1, 2)

    def test_get_celpp_challenge_data_dir_name_path_not_found(self):
        self.assertEqual(1, 2)

    def test_get_celpp_challenge_data_dir_name_not_standard_dir(self):
        self.assertEqual(1, 2)

    def test_get_celpp_challenge_data_dir_name(self):
        self.assertEqual(1, 2)

    def test_create_challenge_dir_unable_to_make_dir(self):
        self.assertEqual(1, 2)

    def test_create_challenge_dir(self):
        self.assertEqual(1, 2)

    def test_create_readme_unable_to_write_file(self):
        self.assertEqual(1, 2)

    def test_create_readme_version_unset(self):
        self.assertEqual(1, 2)

    def test_create_readme_version_no_blastnfilter_summary(self):
        self.assertEqual(1, 2)

    def test_create_readme_version_empty_summary(self):
        self.assertEqual(1, 2)

    def test_create_readme(self):
        self.assertEqual(1, 2)

    def test_copy_over_tsv_files_no_dataimport(self):
        self.assertEqual(1, 2)

    def test_copy_over_tsv_files_no_missing_various_tsv_files(self):
        self.assertEqual(1, 2)

    def test_tar_challenge_dir_unable_to_write_tarfile(self):
        self.assertEqual(1, 2)

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
