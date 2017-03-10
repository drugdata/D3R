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

    def test_get_numdocked_and_average(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test None for Pickle file
            num, avg = post_evaluation.get_numdocked_and_average(None,
                                                                 ctype=None)
            self.assertEqual(num, -1)
            self.assertEqual(avg, -1)

            apickle = os.path.join(temp_dir, 'apickle.pickle')
            data = {'5ka3': {'SMCSS': 1.0, 'hiResHolo': 2.0, 'hiTanimoto': 'foo',
                             post_evaluation.LMCSS: 4.0},
                    '5uud': {'hiResHolo': 5.0, 'SMCSS': 6.0, 'hiResApo': 7.0,
                             'hiTanimoto': 8.0, post_evaluation.LMCSS: 9.0}}
            f = open(apickle, 'w')
            pickle.dump(data, f)
            f.flush()
            f.close()

            # test None for ctype
            num, avg = post_evaluation.get_numdocked_and_average(apickle,
                                                                 ctype=None)
            self.assertEqual(num, -1)
            self.assertEqual(avg, -1)

            # test where no matches found
            num, avg = post_evaluation.get_numdocked_and_average(apickle,
                                                                 ctype='foo')
            self.assertEqual(num, 0)
            self.assertEqual(avg, -1)

            # test where 1 dock is found
            num, avg = post_evaluation.\
                get_numdocked_and_average(apickle,
                                          ctype='hiResApo')
            self.assertEqual(num, 1)
            self.assertEqual(avg, 7.0)

            # test where 2 docked is found
            num, avg = post_evaluation.get_numdocked_and_average(apickle,
                                                                 ctype='SMCSS')
            self.assertEqual(num, 2)
            self.assertEqual(avg, 3.5)

            # test default candidate type
            num, avg = post_evaluation.get_numdocked_and_average(apickle)
            self.assertEqual(num, 2)
            self.assertEqual(avg, 6.5)

            # test where score is not a number
            num, avg = post_evaluation.\
                get_numdocked_and_average(apickle,
                                          ctype='hiTanimoto')
            self.assertEqual(num, -1)
            self.assertEqual(avg, -1)
        finally:
            shutil.rmtree(temp_dir)

    def test_main_success_cause_no_evaluations_to_check(self):
        temp_dir = tempfile.mkdtemp()

        try:
            theargs = ['post_evaluation.py', temp_dir]

            self.assertEqual(post_evaluation.main(theargs), 0)

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
