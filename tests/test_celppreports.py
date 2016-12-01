#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_celppreports
----------------------------------

Tests for `celppreports` module.
"""

import unittest
import tempfile
import logging
import os
import os.path
import shutil

from d3r import celppreports
from d3r.celpp.task import D3RParameters
from d3r.celpp.blastnfilter import BlastNFilterTask


class TestCelppReports(unittest.TestCase):

    def setUp(self):
        pass

    def test_setup_logging(self):
        theargs = D3RParameters()

        theargs.loglevel = 'DEBUG'
        celppreports._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.blastnfilter')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.dataimport')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.glide')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.makeblastdb')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.proteinligprep')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.evaluation')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.task')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.util')
                         .getEffectiveLevel(), logging.DEBUG)
        self.assertEqual(logging.getLogger('d3r.celpp.filetransfer')
                         .getEffectiveLevel(), logging.DEBUG)

        theargs.loglevel = 'INFO'
        celppreports._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.blastnfilter')
                         .getEffectiveLevel(), logging.INFO)

        theargs.loglevel = 'WARNING'
        celppreports._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.blastnfilter')
                         .getEffectiveLevel(), logging.WARNING)

        theargs.loglevel = 'ERROR'
        celppreports._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.blastnfilter')
                         .getEffectiveLevel(), logging.ERROR)

        theargs.loglevel = 'CRITICAL'
        celppreports._setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.blastnfilter')
                         .getEffectiveLevel(), logging.CRITICAL)

    def test_parse_arguments(self):
        theargs = ['--outdir', '/outdir', '--log', 'INFO', 'celppdir']
        result = celppreports._parse_arguments('hi', theargs)
        self.assertEqual(result.outdir, '/outdir')
        self.assertEqual(result.loglevel, 'INFO')
        self.assertEqual(result.celppdir, 'celppdir')

    def test_generate_reports_celppdir_does_not_exist(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()
            theargs.celppdir = os.path.join(temp_dir, 'doesnotexist')
            try:
                celppreports.generate_reports(theargs)
                self.assertEqual('Expected exception')
            except Exception:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_generate_reports_outdir_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()

            theargs.celppdir = temp_dir
            try:
                celppreports.generate_reports(theargs)
                self.assertEqual('Expected exception')
            except Exception:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_generate_reports_outdir_set_to_none(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()

            theargs.celppdir = temp_dir
            theargs.outdir = None
            try:
                celppreports.generate_reports(theargs)
                self.assertEqual('Expected exception')
            except Exception:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_generate_reports_outdir_needs_to_be_created_no_celpp_data(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()

            theargs.celppdir = temp_dir
            outdir = os.path.join(temp_dir, 'outdir')
            theargs.outdir = outdir

            celppreports.generate_reports(theargs)
            self.assertEqual(os.path.isdir(outdir), True)
            csv_file = os.path.join(outdir, 'blastnfilter.summary.csv')
            self.assertEqual(os.path.isfile(csv_file), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_generate_reports_with_one_week_in_one_year(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = D3RParameters()

            theargs.celppdir = temp_dir
            outdir = os.path.join(temp_dir, 'outdir')
            theargs.outdir = outdir

            yeardir = os.path.join(temp_dir, '2016')
            os.mkdir(yeardir)
            weekdir = os.path.join(yeardir, 'dataset.week.10')
            os.mkdir(weekdir)
            blast = BlastNFilterTask(temp_dir, theargs)
            blastdir = os.path.join(weekdir, blast.get_dir_name())
            os.mkdir(blastdir)

            # create dummary summary.txt file

            f = open(os.path.join(blastdir,
                                  BlastNFilterTask.SUMMARY_TXT), 'w')
            f.write('INPUT SUMMARY\n')
            f.write('  entries:                             221\n')
            f.write('  complexes:                           178\n')
            f.write('  dockable complexes:                   95\n')
            f.write('  monomers:                            145\n')
            f.write('  dockable monomers:                    71\n')
            f.write('  multimers:                            76\n')
            f.write('  dockable multimers:                   24\n\n')

            f.write('FILTERING CRITERIA\n')
            f.write('  No. of query sequences           <=    1\n')
            f.write('  No. of dockable ligands           =    1\n')
            f.write('  Percent identity                 >=    0.95\n')
            f.write('  Percent Coverage                 >=    0.9\n')
            f.write('  No. of hit sequences             <=    4\n')
            f.write('  Structure determination method:        '
                    '  x-ray diffraction\n\n')

            f.write('OUTPUT SUMMARY\n')
            f.write('  Targets found:                        67\n')
            f.write('  Target: 5fz7|Sequences: 1|Hits: 94|Candidates: 17|'
                    'Elected:4|PDBids: 5fz7,5fyz,5a3p,5a1f\n')

            f.flush()
            f.close()

            celppreports.generate_reports(theargs)
            self.assertEqual(os.path.isdir(outdir), True)
            csv_file = os.path.join(outdir, 'blastnfilter.summary.csv')
            self.assertEqual(os.path.isfile(csv_file), True)

            # check csv file
            f = open(csv_file, 'r')
            header = f.readline()
            self.assertEqual(header, 'Week #, Year, Complexes, Dockable '
                                     'complexes, Dockable monomers, '
                                     'Targets Found\n')
            data = f.readline()
            self.assertEqual(data, '10,2016,178,95,71,67\n')

        finally:
            shutil.rmtree(temp_dir)

    def test_main_where_generate_reports_raises_error(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = ['prog', temp_dir]
            self.assertEqual(celppreports.main(theargs), 2)
        finally:
            shutil.rmtree(temp_dir)

    def test_main_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            theargs = ['prog', '--outdir', os.path.join(temp_dir, 'foo'),
                       temp_dir]
            self.assertEqual(celppreports.main(theargs), 0)
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
