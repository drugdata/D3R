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

        med_form = '%-15.3f'
        form = '%-15.3f, '
        # pass in list with None and 1 value
        res = evaluate.calculate_average_min_max_median([None, 5])
        avg = form % 5
        min = form % 5
        max = form % 5
        med = med_form % 5
        self.assertEqual(res, (avg, min, max, med))

        # pass in list with 3 values
        res = evaluate.calculate_average_min_max_median([3, 5, 9])
        avg = form % 5.666666
        min = form % 3
        max = form % 9
        med = med_form % 5
        self.assertEqual(res, (avg, min, max, med))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
