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

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = GlideTask(temp_dir, params)
            # try with no dir
            self.assertEqual(task.get_uploadable_files(), [])

            # try with empty dir
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # try with final log
            final_log = os.path.join(task.get_dir(),
                                     GlideTask.FINAL_LOG)
            open(final_log, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(final_log)

            # try with empty pbdid dir
            pbdid = os.path.join(task.get_dir(), '6fff')
            os.mkdir(pbdid)
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            flist.index(final_log)

            # try with pbdid.txt
            pbdidtxt = os.path.join(task.get_dir(), '6fff.txt')
            open(pbdidtxt, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(pbdidtxt)

            # try with LMCSS/LMCSS_dock_pv.maegz
            LMCSSd = os.path.join(pbdid, 'LMCSS')
            os.mkdir(LMCSSd)
            LMCSS = os.path.join(LMCSSd, 'LMCSS_dock_pv.maegz')
            open(LMCSS, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(LMCSS)

            # try with SMCSS/SMCSS_dock_pv.maegz
            SMCSSd = os.path.join(pbdid, 'SMCSS')
            os.mkdir(SMCSSd)
            SMCSS = os.path.join(SMCSSd, 'SMCSS_dock_pv.maegz')
            open(SMCSS, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 4)
            flist.index(SMCSS)

            # try with hiResApo/hiResApo_dock_pv.maegz
            hiResApod = os.path.join(pbdid, 'hiResApo')
            os.mkdir(hiResApod)
            hiResApo = os.path.join(hiResApod, 'hiResApo_dock_pv.maegz')
            open(hiResApo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 5)
            flist.index(hiResApo)

            # try with hiResHolo/hiResHolo_dock_pv.maegz
            hiResHolod = os.path.join(pbdid, 'hiResHolo')
            os.mkdir(hiResHolod)
            hiResHolo = os.path.join(hiResHolod, 'hiResHolo_dock_pv.maegz')
            open(hiResHolo, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 6)
            flist.index(hiResHolo)

            # try with stderr/out files
            errfile = os.path.join(task.get_dir(), 'glidedocking.py.stderr')
            open(errfile, 'a').close()
            outfile = os.path.join(task.get_dir(), 'glidedocking.py.stdout')
            open(outfile, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 8)
            flist.index(errfile)
            flist.index(outfile)
            flist.index(final_log)
            flist.index(LMCSS)
            flist.index(SMCSS)
            flist.index(hiResApo)
            flist.index(hiResHolo)
            flist.index(pbdidtxt)

        finally:
            shutil.rmtree(temp_dir)

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
