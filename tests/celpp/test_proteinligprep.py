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
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.proteinligprep import ProteinLigPrepTask


class TestProteinLigPrepTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = ProteinLigPrepTask(temp_dir, params)
            # try with no dir
            self.assertEqual(task.get_uploadable_files(), [])

            # try with empty dir
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # try with final log
            final_log = os.path.join(task.get_dir(),
                                     ProteinLigPrepTask.FINAL_LOG)
            open(final_log, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(final_log)

            # try with pbid folder that is empty
            pbdid = os.path.join(task.get_dir(), '4abc')
            os.mkdir(pbdid)
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)

            # try with pbdid folder with ligand.mae
            ligand = os.path.join(pbdid, 'ligand.mae')
            open(ligand, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(ligand)

            # try with pbdid folder with largest.maegz
            largest = os.path.join(pbdid, 'largest.maegz')
            open(largest, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(largest)

            # try with pbdid folder with smallest.maegz
            smallest = os.path.join(pbdid, 'smallest.maegz')
            open(smallest, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 4)
            flist.index(smallest)

            # try with pbdid folder with apo.maegz
            apo = os.path.join(pbdid, 'apo.maegz')
            open(apo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 5)
            flist.index(apo)

            # try with pbdid folder with holo.maegz
            holo = os.path.join(pbdid, 'holo.maegz')
            open(holo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 6)
            flist.index(holo)

            # add error out files and try with second pbdid folder
            # with ligand.mae
            errfile = os.path.join(task.get_dir(), 'proteinligprep.py.stderr')
            open(errfile, 'a').close()
            outfile = os.path.join(task.get_dir(), 'proteinligprep.py.stdout')
            open(outfile, 'a').close()

            pbdtwo = os.path.join(task.get_dir(), '3zaa')
            os.mkdir(pbdtwo)
            ligandtwo = os.path.join(pbdtwo, 'ligand.mae')
            open(ligandtwo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 9)
            flist.index(errfile)
            flist.index(outfile)
            flist.index(ligand)
            flist.index(ligandtwo)
            flist.index(apo)
            flist.index(largest)
            flist.index(final_log)
            flist.index(holo)
            flist.index(smallest)

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no blast task found so it cannot run
            params = D3RParameters()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             'challengedata task has notfound status')

            # challenge filter running
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.START_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             'challengedata task has start status')

            # blastnfilter failed
            error_file = os.path.join(chall.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             'challengedata task has error status')

            # blastnfilter success
            os.remove(error_file)
            open(os.path.join(chall.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), True)
            self.assertEqual(proteinligprep.get_error(), None)

            # proteinligprep task exists already
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             proteinligprep.get_dir_name() +
                             ' already exists and status is unknown')

            # proteinlibprep already complete
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            open(os.path.join(proteinligprep.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # return immediately cause can_run is false
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(),
                             'challengedata task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_proteinligprep_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(),
                             'proteinligprep not set')
            # test files get created
            self.assertEqual(os.path.isdir(proteinligprep.get_dir()),
                             True)
            errfile = os.path.join(proteinligprep.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_pdbdb_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.proteinligprep = 'false'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(),
                             'pdbdb not set')
            # test files get created
            self.assertEqual(os.path.isdir(proteinligprep.get_dir()),
                             True)
            errfile = os.path.join(proteinligprep.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_proteinligprep_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.proteinligprep = 'false'
            params.pdbdb = '/foo'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)

            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(proteinligprep.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(proteinligprep.get_dir(),
                                  'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(proteinligprep.get_dir(),
                                  'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_proteinligprep_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.proteinligprep = '/bin/doesnotexist'
            params.pdbdb = '/foo'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)

            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --candidatedir ' +
                             chall.get_dir() + ' --pdbdb ' +
                             '/foo --outdir ' +
                             proteinligprep.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(proteinligprep.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.proteinligprep = 'true'
            params.pdbdb = '/foo'
            chall = ChallengeDataTask(temp_dir, params)
            chall.create_dir()
            open(os.path.join(chall.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)

            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(), None)
            # test files get created
            errfile = os.path.join(proteinligprep.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(proteinligprep.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(proteinligprep.get_dir(),
                                  'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(proteinligprep.get_dir(),
                                  'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
