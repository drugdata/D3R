__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_glide
--------------------------------

Tests for `glide` module.
"""

import shutil
import os
from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.glide import GlideTask
from d3r.celpp.proteinligprep import ProteinLigPrepTask


class TestGlideTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_can_run(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # no proteinligprep task found so it cannot run
            params = D3RParameters()
            glide = GlideTask(temp_dir, params)
            self.assertEqual(glide.can_run(), False)
            self.assertEqual(glide.get_error(),
                             'proteinligprep task has notfound status')

            # proteinligprep filter running
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.START_FILE),
                 'a').close()
            glide = GlideTask(temp_dir, params)
            self.assertEqual(glide.can_run(), False)
            self.assertEqual(glide.get_error(),
                             'proteinligprep task has start status')

            # proteinligprep failed
            error_file = os.path.join(proteinligprep.get_dir(),
                                      D3RTask.ERROR_FILE)
            open(error_file, 'a').close()
            glide = GlideTask(temp_dir, params)
            self.assertEqual(glide.can_run(), False)
            self.assertEqual(glide.get_error(),
                             'proteinligprep task has error status')

            # proteinligprep success
            os.remove(error_file)
            open(os.path.join(proteinligprep.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            glide = GlideTask(temp_dir, params)
            self.assertEqual(glide.can_run(), True)
            self.assertEqual(glide.get_error(), None)

            # glide task exists already
            glide = GlideTask(temp_dir, params)
            glide.create_dir()
            self.assertEqual(glide.can_run(), False)
            self.assertEqual(glide.get_error(),
                             glide.get_dir_name() +
                             ' already exists and status is unknown')

            # glide already complete
            glide = GlideTask(temp_dir, params)
            open(os.path.join(glide.get_dir(),
                              D3RTask.COMPLETE_FILE), 'a').close()
            self.assertEqual(glide.can_run(), False)
            self.assertEqual(glide.get_error(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            # return immediately cause can_run is false
            glide = GlideTask(temp_dir, params)
            glide.run()
            self.assertEqual(glide.get_error(),
                             'proteinligprep task has notfound status')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_glide_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            glide = GlideTask(temp_dir, params)
            glide.run()
            self.assertEqual(glide.get_error(),
                             'glide not set')
            # test files get created
            self.assertEqual(os.path.isdir(glide.get_dir()),
                             True)
            errfile = os.path.join(glide.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_glide_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.glide = 'false'
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            glide = GlideTask(temp_dir, params)

            glide.run()
            self.assertEqual(glide.get_error(),
                             'Non zero exit code: 1 received. Standard out: ' +
                             ' Standard error: ')
            # test file gets created
            errfile = os.path.join(glide.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)

            stderr = os.path.join(glide.get_dir(),
                                   'false.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(glide.get_dir(),
                                   'false.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_fails_cause_glide_is_not_found(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.glide = '/bin/doesnotexist'
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()
            glide = GlideTask(temp_dir, params)

            glide.run()
            self.assertEqual(glide.get_error(),
                             'Caught Exception trying to run ' +
                             '/bin/doesnotexist --structuredir ' +
                             proteinligprep.get_dir() + ' --outdir ' +
                             glide.get_dir() +
                             ' : [Errno 2] No such file or directory')

            # test files get created
            errfile = os.path.join(glide.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.glide = 'true'
            proteinligprep = ProteinLigPrepTask(temp_dir, params)
            proteinligprep.create_dir()
            open(os.path.join(proteinligprep.get_dir(),
                              D3RTask.COMPLETE_FILE),
                 'a').close()
            glide = GlideTask(temp_dir, params)

            glide.run()
            self.assertEqual(glide.get_error(), None)
            # test files get created
            errfile = os.path.join(glide.get_dir(),
                                  D3RTask.ERROR_FILE)
            self.assertEqual(os.path.isfile(errfile), False)

            compfile = os.path.join(glide.get_dir(),
                                    D3RTask.COMPLETE_FILE)
            self.assertEqual(os.path.isfile(compfile), True)
            stderr = os.path.join(glide.get_dir(),
                                   'true.stderr')
            self.assertEqual(os.path.isfile(stderr), True)
            stdout = os.path.join(glide.get_dir(),
                                   'true.stdout')
            self.assertEqual(os.path.isfile(stdout), True)
        finally:
            shutil.rmtree(temp_dir)
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()