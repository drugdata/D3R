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
from d3r.celpp.proteinligprep import ProteinLigPrepTask


class TestProteinLigPrepTask(unittest.TestCase):
    def setUp(self):
        pass


    def test_get_uploadable_files(self):
        self.assertEqual('need to ', 'test this')

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no blast task found so it cannot run
            params = D3RParameters()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             'blastnfilter task has notfound status')

            # blastn filter running
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.START_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             'blastnfilter task has start status')

            # blastnfilter failed
            error_file = os.path.join(blastnfilter.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            self.assertEqual(proteinligprep.can_run(), False)
            self.assertEqual(proteinligprep.get_error(),
                             'blastnfilter task has error status')

            # blastnfilter success
            os.remove(error_file)
            open(os.path.join(blastnfilter.get_dir(),
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
                             'blastnfilter task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_proteinligprep_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
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
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
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
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
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
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)

            proteinligprep.run()
            self.assertEqual(proteinligprep.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --candidatedir ' +
                             blastnfilter.get_dir() + ' --pdbdb ' +
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
            blastnfilter = BlastNFilterTask(temp_dir, params)
            blastnfilter.create_dir()
            open(os.path.join(blastnfilter.get_dir(), D3RTask.COMPLETE_FILE),
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
