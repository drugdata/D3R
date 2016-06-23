__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_vina
--------------------------------

Tests for `vina` module.
"""

import shutil
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.vina import AutoDockVinaTask
from d3r.celpp.chimeraprep import ChimeraProteinLigPrepTask


class TestAutoDockVinaTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = AutoDockVinaTask(temp_dir, params)
            # try with no dir
            self.assertEqual(task.get_uploadable_files(), [])

            # try with empty dir
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # try with stderr/out files
            errfile = os.path.join(task.get_dir(), 'vinadocking.py.stderr')
            open(errfile, 'a').close()
            outfile = os.path.join(task.get_dir(), 'vinadocking.py.stdout')
            open(outfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(errfile)
            flist.index(outfile)

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no proteinligprep task found so it cannot run
            params = D3RParameters()
            vina = AutoDockVinaTask(temp_dir, params)
            self.assertEqual(vina.can_run(), False)
            self.assertEqual(vina.get_error(),
                             'chimeraprep task has notfound status')

            # proteinligprep filter running
            proteinligprep = ChimeraProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.START_FILE),
                 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)
            self.assertEqual(vina.can_run(), False)
            self.assertEqual(vina.get_error(),
                             'chimeraprep task has start status')

            # proteinligprep failed
            error_file = os.path.join(proteinligprep.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)
            self.assertEqual(vina.can_run(), False)
            self.assertEqual(vina.get_error(),
                             'chimeraprep task has error status')

            # proteinligprep success
            os.remove(error_file)
            open(os.path.join(proteinligprep.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)
            self.assertEqual(vina.can_run(), True)
            self.assertEqual(vina.get_error(), None)

            # vina task exists already
            vina = AutoDockVinaTask(temp_dir, params)
            vina.create_dir()
            self.assertEqual(vina.can_run(), False)
            self.assertEqual(vina.get_error(),
                             vina.get_dir_name() +
                             ' already exists and status is unknown')

            # vina already complete
            vina = AutoDockVinaTask(temp_dir, params)
            open(os.path.join(vina.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(vina.can_run(), False)
            self.assertEqual(vina.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # return immediately cause can_run is false
            vina = AutoDockVinaTask(temp_dir, params)
            vina.run()
            self.assertEqual(vina.get_error(),
                             'chimeraprep task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_vina_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            proteinligprep = ChimeraProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)
            vina.run()
            self.assertEqual(vina.get_error(),
                             'vina not set')
            # test files get created
            self.assertEqual(os.path.isdir(vina.get_dir()),
                             True)
            errfile = os.path.join(vina.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_vina_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.vina = 'false'
            proteinligprep = ChimeraProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)

            vina.run()
            self.assertEqual(vina.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(vina.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(vina.get_dir(),
                                  'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(vina.get_dir(),
                                  'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_vina_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.vina = '/bin/doesnotexist'
            proteinligprep = ChimeraProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)

            vina.run()
            self.assertEqual(vina.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --structuredir ' +
                             proteinligprep.get_dir() + ' --outdir ' +
                             vina.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(vina.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.vina = 'true'
            proteinligprep = ChimeraProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(),
                              D3RTask.COMPLETE_FILE),
                 'a').close()
            vina = AutoDockVinaTask(temp_dir, params)

            vina.run()
            self.assertEqual(vina.get_error(), None)
            # test files get created
            errfile = os.path.join(vina.get_dir(),
                                   D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(vina.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(vina.get_dir(),
                                  'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(vina.get_dir(),
                                  'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
