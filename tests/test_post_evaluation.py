#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_post_evaluation
----------------------------------

Tests for `post_evaluation` module.
"""

import unittest
import tempfile
import os
import os.path
import shutil
import pickle

from d3r import post_evaluation
from d3r.celpp.task import D3RParameters


class TestPostEvaluation(unittest.TestCase):
    """Tests post_evaluation commandline script
    """
    def setUp(self):
        pass

    def test_check_case_number(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # log file is none
            self.assertEqual(post_evaluation.check_case_number(None,
                                                               'hi'), -1)

            # non existant file
            noexist = os.path.join(temp_dir, 'doesnotexist.txt')
            self.assertEqual(post_evaluation.check_case_number(noexist,
                                                               'hi'), -1)

            # passing directory which should fail
            self.assertEqual(post_evaluation.check_case_number(temp_dir,
                                                               'hi'), -1)

            # empty file
            efile = os.path.join(temp_dir, 'emptyfile.txt')
            open(efile, 'a').close()
            self.assertEqual(post_evaluation.check_case_number(efile,
                                                               'hi'), 0)

            # phrase is None
            self.assertEqual(post_evaluation.check_case_number(efile,
                                                               None), -1)
            # file with couple matches
            tfile = os.path.join(temp_dir, 'hi.txt')
            f = open(tfile, 'w')
            f.write('Hi\n bye some phrase\n some\n'
                    ' phrase\ some phrase sdfsdf some phrase\n')
            f.flush()
            f.close()
            self.assertEqual(post_evaluation.check_case_number(tfile,
                                                               'some phrase'),
                             2)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_dock_scores_as_list(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test passing dir as pickle file
            res = post_evaluation.get_dock_scores_as_list(temp_dir)
            self.assertEqual(res, [])

            # test None for Pickle file
            res = post_evaluation.get_dock_scores_as_list(None,
                                                          ctype=None)
            self.assertEqual(res, [])

            smcss = post_evaluation.SMCSS
            hi_tanimoto = post_evaluation.HI_TANIMOTO
            hi_resapo = post_evaluation.HI_RESAPO
            hi_resholo = post_evaluation.HI_RESHOLO
            apickle = os.path.join(temp_dir, 'apickle.pickle')
            data = {'5ka3': {smcss: 1.0, hi_resholo: 2.0, hi_tanimoto: 'foo',
                             post_evaluation.LMCSS: 4.0},
                    '5uud': {hi_resholo: 5.0, smcss: 6.0, hi_resapo: 7.0,
                             hi_tanimoto: 8.0, post_evaluation.LMCSS: 9.0}}
            f = open(apickle, 'w')
            pickle.dump(data, f)
            f.flush()
            f.close()

            # test None for ctype
            res = post_evaluation.get_dock_scores_as_list(apickle,
                                                               ctype=None)
            self.assertEqual(res, [])

            # test where no matches found
            res = post_evaluation.get_dock_scores_as_list(apickle,
                                                          ctype='foo')
            self.assertEqual(res, [])

            # test where 1 dock is found
            res = post_evaluation.get_dock_scores_as_list(apickle,
                                                          ctype=hi_resapo)
            self.assertEqual(res, [7.0])

            # test where 2 docked is found
            res = post_evaluation.get_dock_scores_as_list(apickle,
                                                          ctype=smcss)
            self.assertTrue(6.0 in res)
            self.assertTrue(1.0 in res)
            self.assertEqual(len(res), 2)

            # test default candidate type
            res = post_evaluation.get_dock_scores_as_list(apickle)
            self.assertTrue(4.0 in res)
            self.assertTrue(9.0 in res)
            self.assertEqual(len(res), 2)

            # test where score is not a number
            res = post_evaluation.get_dock_scores_as_list(apickle,
                                                          ctype=hi_tanimoto)
            self.assertTrue('foo' in res)
            self.assertTrue(8.0 in res)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_list_of_stats(self):
        res = post_evaluation.get_list_of_stats(None)
        self.assertEqual(res, (-1, -1, -1, -1))

        res = post_evaluation.get_list_of_stats([1.0])
        self.assertEqual(res, (1, 1.0, 1.0, 1.0))

        res = post_evaluation.get_list_of_stats([3.0, 1.0])
        self.assertEqual(res, (2, 1.0, 3.0, 2.0))

        res = post_evaluation.get_list_of_stats([1.0, 'foo'])
        self.assertEqual(res, (-1, -1, -1, -1))

        res = post_evaluation.get_list_of_stats([3.0, 1.0, 2.0, 14.0])
        self.assertEqual(res, (4, 1.0, 14.0, 5.0))

    def test_get_pickle_paths(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # none passed in
            p_list, count = post_evaluation._get_pickle_paths(None)
            self.assertEqual(p_list, [])
            self.assertEqual(count, 0)

            # emptylist
            p_list, count = post_evaluation._get_pickle_paths([])
            self.assertEqual(p_list, [])
            self.assertEqual(count, 0)

            # one path no pickle
            p_list, count = post_evaluation._get_pickle_paths([temp_dir])
            self.assertEqual(p_list, [])
            self.assertEqual(count, 1)

            # one path with pickle
            pfile = os.path.join(temp_dir, post_evaluation.RMSD_PICKLE)
            open(pfile, 'a').close()
            p_list, count = post_evaluation._get_pickle_paths([temp_dir])
            self.assertEqual(len(p_list), 1)
            self.assertTrue(pfile in p_list)
            self.assertEqual(count, 0)

            # two paths and only one with pickle
            subdir = os.path.join(temp_dir, 'foo')
            os.makedirs(subdir, mode=0o775)
            p_list, count = post_evaluation._get_pickle_paths([temp_dir,
                                                               subdir])
            self.assertEqual(len(p_list), 1)
            self.assertTrue(pfile in p_list)
            self.assertEqual(count, 1)

            # two paths with pickles
            p2file = os.path.join(subdir, post_evaluation.RMSD_PICKLE)
            open(p2file, 'a').close()

            p_list, count = post_evaluation._get_pickle_paths([subdir,
                                                               temp_dir])
            self.assertEqual(len(p_list), 2)
            self.assertTrue(pfile in p_list)
            self.assertTrue(p2file in p_list)
            self.assertEqual(count, 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_submission_name_from_pickle_path(self):

        somepath = os.path.join('data', 'celpp', '2017', 'dataset.week.10',
                                'stage.7.foo_dock.extsubmission.evaluation',
                                post_evaluation.RMSD_PICKLE)
        prefix = 'stage.7.'
        suffix = '.extsubmission.evaluation$|.evaluation$'
        name = post_evaluation._get_submission_name_from_pickle_path(somepath,
                                                                     prefix,
                                                                     suffix)
        self.assertEqual(name, 'foo_dock')

        # non external submission
        somepath = os.path.join('data', 'celpp', '2017', 'dataset.week.10',
                                'stage.7.foo_dock.evaluation',
                                post_evaluation.RMSD_PICKLE)

        self.assertEqual(name, 'foo_dock')

    def test_generate_overall_csv_valid_single_evaldataset(self):
        temp_dir = tempfile.mkdtemp()
        try:
            result_dir = os.path.join(temp_dir, 'result')
            os.makedirs(result_dir, mode=0o755)

            chall_dir = os.path.join(temp_dir, 'challdir')
            os.makedirs(chall_dir, mode=0o755)

            # write out the final.log file
            finalfile = os.path.join(chall_dir, post_evaluation.CHALL_FINAL_LOG)
            f = open(finalfile, 'w')
            f.write("""03/04/17 03:32:54: ============Start to work in this query protein: 5uub============
03/04/17 03:32:55: Warning: Found multiple ligand for this protein:LMCSS-5uub_5un3-XYP.pdb
03/04/17 03:32:55: Succsessfully generate this protein:LMCSS-5uub_5un3-XYP.pdb
03/04/17 03:32:55: Ligand center for this case:LMCSS-5uub_5un3-XYP-lig.pdb is   38.448,   73.584,   58.680
03/04/17 03:33:00: Succsessfully generate this protein:SMCSS-5uub_5fsj-OXY.pdb
03/04/17 03:33:12: ============Start to work in this query protein: 5ur3============
03/04/17 03:33:12: Succsessfully generate this protein:LMCSS-5ur3_3njq-NJQ.pdb
03/04/17 03:33:12: Ligand center for this case:LMCSS-5ur3_3njq-NJQ-lig.pdb is  -27.326,  -19.578,  -43.347
03/04/17 03:33:19: Succsessfully generate this protein:hiResHolo-5ur3_4p3h-25G.pdb
03/04/17 03:33:25: Succsessfully generate this protein:LMCSS-5uud_1qf2-TI3.pdb
03/04/17 03:33:25: Ligand center for this case:LMCSS-5uud_1qf2-TI3-lig.pdb is   38.308,   38.885,   -3.382
03/04/17 03:33:40: Succsessfully generate this protein:hiTanimoto-5uud_1gxw-VAL.pdb""")
            f.flush()
            f.close()

            task_dir = os.path.join(temp_dir,
                                    'stage.7.foo_dock.extsubmission.evaluation')
            os.makedirs(task_dir, mode=0o755)

            # write out the pickle dir
            smcss = post_evaluation.SMCSS
            hi_tanimoto = post_evaluation.HI_TANIMOTO
            hi_resapo = post_evaluation.HI_RESAPO
            hi_resholo = post_evaluation.HI_RESHOLO
            apickle = os.path.join(task_dir, post_evaluation.RMSD_PICKLE)
            data = {'5ka3': {smcss: 1.0, hi_resholo: 2.0, hi_tanimoto: 'foo',
                             post_evaluation.LMCSS: 4.0},
                    '5uud': {hi_resholo: 5.0, smcss: 6.0, hi_resapo: 7.0,
                             hi_tanimoto: 8.0, post_evaluation.LMCSS: 9.0}}
            f = open(apickle, 'w')
            pickle.dump(data, f)
            f.flush()
            f.close()
            post_evaluation.generate_overall_csv([task_dir], chall_dir,
                                                 result_dir,
                                                 post_evaluation.LMCSS)

            summary_file = os.path.join(result_dir,
                                        post_evaluation.SUMMARY_TXT)
            self.assertTrue(os.path.isfile(summary_file))
            f = open(summary_file, 'r')
            data = f.read()
            f.close()

            self.assertTrue(post_evaluation.LMCSS + ' (3 dockable)' in data)
            self.assertTrue(('stage.7.foo_dock.' +
                            'extsubmission.evaluation') in data)
            self.assertTrue('2 ( 67%)    4.000       9.000' in data)

            csv_file = os.path.join(result_dir,
                                    post_evaluation.OVERALL_RMSD_PREFIX +
                                    post_evaluation.LMCSS +
                                    post_evaluation.CSV_SUFFIX)
            self.assertTrue(os.path.isfile(csv_file))
            f = open(csv_file, 'r')
            data = f.read()
            f.close()
            self.assertTrue('SubmissionID for LMCSS,# Docked,' in data)
            self.assertTrue('stage.7.foo_dock.extsubmission' in data)
            self.assertTrue('submission.evaluation,2,3,4.000,9.000,6.500' in data)

            # test prefix suffix removal
            os.remove(summary_file)
            post_evaluation.generate_overall_csv([task_dir], chall_dir,
                                                 result_dir,
                                                 post_evaluation.LMCSS,
                                                 eval_stage_prefix='stage.7.',
                                                 eval_suffix='.extsubmission.'
                                                             'evaluation')

            f = open(summary_file, 'r')
            data = f.read()
            f.close()
            self.assertTrue('foo_dock ' in data)
            f = open(csv_file, 'r')
            data = f.read()
            f.close()
            self.assertTrue('foo_dock,2,3,4.000,9.000,6.500' in data)
        finally:
            shutil.rmtree(temp_dir)

    def test_main_success_cause_no_evaluations_to_check(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = ['post_evaluation.py', temp_dir]

            self.assertEqual(post_evaluation.main(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_main_with_eval_task(self):
        temp_dir = tempfile.mkdtemp()

        try:
            result_dir = os.path.join(temp_dir, 'result')
            os.makedirs(result_dir, mode=0o755)

            chall_dir = os.path.join(temp_dir, 'challdir')
            os.makedirs(chall_dir, mode=0o755)

            # write out the final.log file
            finalfile = os.path.join(chall_dir, post_evaluation.CHALL_FINAL_LOG)
            f = open(finalfile, 'w')
            f.write("""03/04/17 03:32:54: ============Start to work in this query protein: 5uub============
            03/04/17 03:32:55: Warning: Found multiple ligand for this protein:LMCSS-5uub_5un3-XYP.pdb
            03/04/17 03:32:55: Succsessfully generate this protein:LMCSS-5uub_5un3-XYP.pdb
            03/04/17 03:32:55: Ligand center for this case:LMCSS-5uub_5un3-XYP-lig.pdb is   38.448,   73.584,   58.680
            03/04/17 03:33:00: Succsessfully generate this protein:SMCSS-5uub_5fsj-OXY.pdb
            03/04/17 03:33:12: ============Start to work in this query protein: 5ur3============
            03/04/17 03:33:12: Succsessfully generate this protein:LMCSS-5ur3_3njq-NJQ.pdb
            03/04/17 03:33:12: Ligand center for this case:LMCSS-5ur3_3njq-NJQ-lig.pdb is  -27.326,  -19.578,  -43.347
            03/04/17 03:33:19: Succsessfully generate this protein:hiResHolo-5ur3_4p3h-25G.pdb
            03/04/17 03:33:25: Ligand center for this case:LMCSS-5uud_1qf2-TI3-lig.pdb is   38.308,   38.885,   -3.382
            03/04/17 03:33:40: Succsessfully generate this protein:hiTanimoto-5uud_1gxw-VAL.pdb""")
            f.flush()
            f.close()

            task_dir = os.path.join(temp_dir,
                                    'stage.7.foo_dock.extsubmission.evaluation')
            os.makedirs(task_dir, mode=0o755)

            # write out the pickle dir
            smcss = post_evaluation.SMCSS
            hi_tanimoto = post_evaluation.HI_TANIMOTO
            hi_resapo = post_evaluation.HI_RESAPO
            hi_resholo = post_evaluation.HI_RESHOLO
            apickle = os.path.join(task_dir, post_evaluation.RMSD_PICKLE)
            data = {'5ka3': {smcss: 1.0, hi_resholo: 2.0, hi_tanimoto: 'foo',
                             post_evaluation.LMCSS: 4.0},
                    '5uud': {hi_resholo: 5.0, smcss: 6.0, hi_resapo: 7.0,
                             hi_tanimoto: 8.0, post_evaluation.LMCSS: 9.0}}
            f = open(apickle, 'w')
            pickle.dump(data, f)
            f.flush()
            f.close()
            theargs = ['post_evaluation.py', result_dir,
                       '--evaluationdir', task_dir,
                       '--challengedir', chall_dir]

            self.assertEqual(post_evaluation.main(theargs), 0)

            summary_file = os.path.join(result_dir,
                                        post_evaluation.SUMMARY_TXT)
            self.assertTrue(os.path.isfile(summary_file))
            f = open(summary_file, 'r')
            data = f.read()
            f.close()

            self.assertTrue(post_evaluation.LMCSS + ' (2 dockable)' in data)
            self.assertTrue('foo_dock ' in data)
            self.assertTrue(' 2 (100%)    4.000       9.000' in data)

            csv_file = os.path.join(result_dir,
                                    post_evaluation.OVERALL_RMSD_PREFIX +
                                    post_evaluation.LMCSS +
                                    post_evaluation.CSV_SUFFIX)
            self.assertTrue(os.path.isfile(csv_file))
            f = open(csv_file, 'r')
            data = f.read()
            f.close()
            self.assertTrue('SubmissionID for LMCSS,# Docked,' in data)
            self.assertTrue('foo_dock,2,2,4.000,9.000,6.500' in data)

            # TODO check contents of other CSV files
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
