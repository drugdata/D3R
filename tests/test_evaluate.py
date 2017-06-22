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

from d3r import evaluate
from d3r.evaluate import data_container


class TestEvaluate(unittest.TestCase):
    """Tests post_evaluation commandline script
    """
    def setUp(self):
        pass

    def test_clean_up_list_of_value(self):
        # test on None
        res = evaluate.clean_up_list_of_value(None)
        self.assertEqual(res, [])

        # test on empty list
        res = evaluate.clean_up_list_of_value([])
        self.assertEqual(res, [])

        # test on list with only None values
        res = evaluate.clean_up_list_of_value([None, None])
        self.assertEqual(res, [])

        # test on list mixed list
        res = evaluate.clean_up_list_of_value([None, 1, 2, None])
        self.assertEqual(res, [1, 2])

    def test_calculate_average_min_max_median(self):
        # pass in None
        res = evaluate.calculate_average_min_max_median(None)
        spacefield = "%-15s, " % (" ")
        self.assertEqual(res, (spacefield, spacefield, spacefield, spacefield))

        # pass in empty list
        res = evaluate.calculate_average_min_max_median([])
        self.assertEqual(res, (spacefield, spacefield, spacefield, spacefield))

        form = '%-15.3f, '
        # pass in list with None and 1 value
        res = evaluate.calculate_average_min_max_median([None, 5])
        avg = form % 5
        min = form % 5
        max = form % 5
        med = form % 5
        self.assertEqual(res, (avg, min, max, med))

        # pass in list with 3 values
        res = evaluate.calculate_average_min_max_median([3, 5, 9])
        avg = form % 5.666666
        min = form % 3
        max = form % 9
        med = form % 5
        self.assertEqual(res, (avg, min, max, med))

    def test_data_container_layout_plain_empty(self):
        temp_dir = tempfile.mkdtemp()
        try:
            dc = data_container()
            rmsd = os.path.join(temp_dir, 'rmsd')
            dc.layout_plain(plain_filename=rmsd)

            rmsdcsv = rmsd + '.csv'
            rmsdtxt = rmsd + '.txt'
            self.assertTrue(os.path.isfile(rmsdcsv))
            self.assertTrue(os.path.isfile(rmsdtxt))

            f = open(rmsdcsv, 'r')
            lines = f.readlines()
            f.close()
            self.assertTrue(lines[0].startswith('Target_PDBID'))
            self.assertTrue(lines[2].startswith('Summary Statistics'))
            self.assertTrue(lines[4].startswith('Number_of_cases'))
            self.assertTrue(lines[5].startswith('Average'))
            self.assertTrue(lines[6].startswith('Maximum'))
            self.assertTrue(lines[7].startswith('Minimum'))
            self.assertTrue(lines[8].startswith('Median'))
            self.assertTrue(lines[10].startswith('Individual Results'))

            f = open(rmsdtxt, 'r')
            lines = f.readlines()
            f.close()
            self.assertTrue(lines[0].startswith('Target_PDBID'))
            self.assertTrue(lines[2].startswith('Summary Statistics'))
            self.assertTrue(lines[4].startswith('Number_of_cases'))
            self.assertTrue(lines[5].startswith('Average'))
            self.assertTrue(lines[6].startswith('Maximum'))
            self.assertTrue(lines[7].startswith('Minimum'))
            self.assertTrue(lines[8].startswith('Median'))
            self.assertTrue(lines[10].startswith('Individual Results'))

        finally:
            shutil.rmtree(temp_dir)

    def test_data_container_layout_plain_with_data(self):
        temp_dir = tempfile.mkdtemp()
        try:
            dc = data_container()
            dc.register('1abc', 'LMCSS', 10.5)
            dc.register('1abc', 'SMCSS', 1.0)
            dc.register('1abc', 'hiResApo', 2.0)
            dc.register('1abc', 'hiResHolo', 3.0)
            dc.register('1abc', 'hiTanimoto', 4.0)
            dc.register('1abc', 'LMCSS_ori', 5.0)
            dc.register('1abc', 'LMCSS_ori_dis', 6.0)

            dc.register('1abc', 'LMCSS_dis', 1.0)
            dc.register('2abc', 'LMCSS', 20.5)
            dc.register('3abc', 'SMCSS', 2.0)
            dc.register('4abc', 'SMCSS', 3.0)


            rmsd = os.path.join(temp_dir, 'rmsd')
            dc.layout_plain(plain_filename=rmsd)

            rmsdcsv = rmsd + '.csv'
            rmsdtxt = rmsd + '.txt'
            self.assertTrue(os.path.isfile(rmsdcsv))
            self.assertTrue(os.path.isfile(rmsdtxt))

            f = open(rmsdcsv, 'r')
            lines = f.readlines()
            f.close()
            self.assertTrue(lines[0].startswith('Target_PDBID'))
            self.assertTrue(lines[2].startswith('Summary Statistics'))
            self.assertTrue(lines[4].startswith('Number_of_cases'))
            self.assertTrue('2              , 3              , 1'
                            '              , 1              , 1' in lines[4])
            self.assertTrue(lines[5].startswith('Average'))
            self.assertTrue('15.500         , 2.000          , 2.000'
                            '          , 3.000          , 4.000' in lines[5])
            self.assertTrue(lines[6].startswith('Maximum'))
            self.assertTrue('20.500         , 3.000          , 2.000'
                            '          , 3.000          , 4.000' in lines[6])
            self.assertTrue(lines[7].startswith('Minimum'))
            self.assertTrue('10.500         , 1.000          , 2.000'
                            '          , 3.000          , 4.000' in lines[7])
            self.assertTrue(lines[8].startswith('Median'))
            self.assertTrue('20.500         , 2.000          , 2.000'
                            '          , 3.000          , 4.000' in lines[8])
            self.assertTrue(lines[10].startswith('Individual Results'))
            self.assertTrue('4abc,                              , 3.000' in lines[12])
            self.assertTrue('3abc,                              , 2.000 ' in lines[13])
            self.assertTrue('2abc,               20.500' in lines[14])
            self.assertTrue('1abc,               10.500(1.000 ) , 1.000'
                            '          , 2.000          , 3.000'
                            '          , 4.000          , 5.000 (6.000 )' in lines[15])
            f = open(rmsdtxt, 'r')
            lines = f.readlines()
            f.close()
            self.assertTrue(lines[0].startswith('Target_PDBID'))
            self.assertTrue(lines[2].startswith('Summary Statistics'))
            self.assertTrue(lines[4].startswith('Number_of_cases'))
            self.assertTrue('2                3                1'
                            '                1                1' in lines[4])
            self.assertTrue(lines[5].startswith('Average'))
            self.assertTrue(lines[6].startswith('Maximum'))
            self.assertTrue(lines[7].startswith('Minimum'))
            self.assertTrue(lines[8].startswith('Median'))
            self.assertTrue(lines[10].startswith('Individual Results'))

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
