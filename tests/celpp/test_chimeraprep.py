__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_chimeraprep
--------------------------------

Tests for `chimeraprep` module.
"""

import shutil
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.chimeraprep import ChimeraProteinLigPrepTask


class TestChimeraProteinLigPrepTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            # try with no dir
            self.assertEqual(task.get_uploadable_files(), [])

            # try with empty dir
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # try with final log
            final_log = os.path.join(task.get_dir(),
                                     ChimeraProteinLigPrepTask.FINAL_LOG)
            open(final_log, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(final_log)

            # add error out files and try with second pbdid folder
            # with ligand.mae
            errfile = os.path.join(task.get_dir(),
                                   'chimera_proteinligprep.py.stderr')
            open(errfile, 'a').close()
            outfile = os.path.join(task.get_dir(),
                                   'chimera_proteinligprep.py.stdout')
            open(outfile, 'a').close()

            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(errfile)
            flist.index(outfile)
            flist.index(final_log)

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no blastnfilter task found so it cannot run
            params = D3RParameters()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             'challengedata task has notfound status')

            # blastnfilter  running
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.START_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             'challengedata task has start status')

            # blastnfilter failed
            error_file = os.path.join(chall.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             'challengedata task has error status')

            # blastnfilter success
            os.remove(error_file)
            open(os.path.join(chall.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            self.assertEqual(task.can_run(), True)
            self.assertEqual(task.get_error(), None)

            # chimeraproteinligprep task exists already
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            task.create_dir()
            self.assertEqual(task.can_run(), False)
            self.assertEqual(task.get_error(),
                             task.get_dir_name() +
                             ' already exists and status is unknown')

            # chimeraproteinligprep already complete
            task = ChimeraProteinLigPrepTask(temp_dir, params)
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
            # return immediately cause can_run is false
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(),
                             'challengedata task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_chimeraprep_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(),
                             'chimeraprep not set')
            # test files get created
            self.assertEqual(os.path.isdir(task.get_dir()),
                             True)
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_pdbdb_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.chimeraprep = 'false'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(),
                             'pdbdb not set')
            # test files get created
            self.assertEqual(os.path.isdir(task.get_dir()),
                             True)
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_chimeraprep_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.chimeraprep = 'false'
            params.pdbdb = '/foo'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)

            task.run()
            self.assertEqual(task.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(task.get_dir(),
                                  'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(task.get_dir(),
                                  'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_chimeraprep_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.chimeraprep = '/bin/doesnotexist'
            params.pdbdb = '/foo'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()

            challdir = os.path.join(chall.get_dir(),
                                    chall.get_celpp_challenge_data_dir_name())

            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)

            task.run()
            self.assertEqual(task.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --candidatedir ' +
                             challdir +
                             ' --rdkitpython \'\' ' +
                             '--pdbdb ' +
                             '/foo --outdir ' +
                             task.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_chimeraprep_is_not_found_rdkitpython_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.chimeraprep = '/bin/doesnotexist'
            params.pdbdb = '/foo'
            params.rdkitpython = '/data/miniconda2/bin'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()

            challdir = os.path.join(chall.get_dir(),
                                    chall.get_celpp_challenge_data_dir_name())

            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)

            task.run()
            self.assertEqual(task.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --candidatedir ' +
                             challdir +
                             ' --rdkitpython \'/data/miniconda2/bin\' ' +
                             '--pdbdb ' +
                             '/foo --outdir ' +
                             task.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.chimeraprep = 'true'
            params.pdbdb = '/foo'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            task = ChimeraProteinLigPrepTask(temp_dir, params)

            task.run()
            self.assertEqual(task.get_error(), None)
            # test files get created
            errfile = os.path.join(task.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(task.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(task.get_dir(),
                                  'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(task.get_dir(),
                                  'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
