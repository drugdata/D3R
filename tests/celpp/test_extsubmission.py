#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'

import unittest
import tempfile
import shutil
import os

from d3r.celpp.extsubmission import ExternalDataSubmissionFactory
from d3r.celpp.task import D3RParameters


"""
test_extsubmission
--------------------------------

Tests for `extsubmission` module.
"""


class TestExternalSubmission(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_externaldatasubmissionfactory_noftpconfig(self):
        params = D3RParameters()
        fac = ExternalDataSubmissionFactory('/foo', params)
        self.assertEqual(fac.get_file_transfer(), None)
        self.assertEqual(fac._get_challenge_package_results_file_name('xxx'),
                         'celpp_week0_0_dockedresults_xxx.tar.gz')

    def test_externaldatasubmissionfactory_ftpconfig_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            figfile = os.path.join(temp_dir,'foo')
            f = open(figfile, 'w')
            f.write('host blah.blah.com\n')
            f.write('user bob\n')
            f.write('pass ha\n')
            f.write('path /celppweekly/ha\n')
            f.write('challengepath /celppweekly/ha\n')
            f.flush()
            f.close()
            params = D3RParameters()
            params.ftpconfig = figfile
            fac = ExternalDataSubmissionFactory('/foo', params)
            self.assertTrue(fac.get_file_transfer() != None)
        finally:
            shutil.rmtree(temp_dir)

    def test_externaldatasubmissionfactory_getter_setter(self):
        params = D3RParameters()
        params.ha = 'hi'
        fac = ExternalDataSubmissionFactory('/foo', params)
        self.assertEqual(fac.get_file_transfer(), None)
        self.assertEqual(fac.get_path(), '/foo')
        self.assertEqual(fac.get_args().ha, 'hi')

        fac.set_file_transfer('yo')
        self.assertEqual(fac.get_file_transfer(), 'yo')


if __name__ == '__main__':
    unittest.main()
