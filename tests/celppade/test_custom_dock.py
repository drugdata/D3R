#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_custom_dock
----------------------------------

Tests for `custom_dock` module.
"""

import unittest
import tempfile
import shutil
import os
import time

from d3r.celppade.custom_dock import Dock


class TestDock(unittest.TestCase):

    def setUp(self):
        pass

    def test_dock(self):
        """Test where we run dock() with fake arguments
           to verify we always get False
        """
        d = Dock()
        self.assertEqual(d.dock(1, 2, 3, 4), False)

    def test_ligand_technical_prep(self):
        """Test where we run ligand_technical_prep() with fake argument
           to verify we always get back the argument as an array
        """
        d = Dock()
        self.assertEqual(d.ligand_technical_prep('hi'), ['hi'])

    def test_receptor_technical_prep(self):
        """Test where we run receptor_technical_prep() with fake arguments
           to verify we always get back the 1st argument as an array
        """
        d = Dock()
        self.assertEqual(d.receptor_technical_prep('hi', 'foo'), ['hi'])

    def test_get_pocket_center(self):
        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            # test where file does not exist
            self.assertEqual(d.get_pocket_center(temp_dir), False)
        finally:
            shutil.rmtree(temp_dir)    


        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            # test with valid center.txt file
            pcfile = os.path.join(temp_dir, 'center.txt')
            f = open(pcfile, 'w')
            f.write('1,2')
            f.flush()
            f.close()
            self.assertEqual(d.get_pocket_center(temp_dir), False)
        finally:
            shutil.rmtree(temp_dir)    


        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            # test where file does not exist
            self.assertEqual(d.get_pocket_center(temp_dir), False)
            pcfile = os.path.join(temp_dir, 'center.txt')
            f = open(pcfile, 'w')
            f.write('1,2,3')
            f.flush()
            f.close()
            self.assertEqual(d.get_pocket_center(temp_dir), [1.0, 2.0, 3.0])
        finally:
            shutil.rmtree(temp_dir)    


        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            # test where file does not exist
            self.assertEqual(d.get_pocket_center(temp_dir), False)
            pcfile = os.path.join(temp_dir, 'center.txt')
            f = open(pcfile, 'w')
            f.write('a,b,c')
            f.flush()
            f.close()
            self.assertEqual(d.get_pocket_center(temp_dir), False)
        finally:
            shutil.rmtree(temp_dir)    

    def test_get_sci_prepped_lig(self):
        # test where file does not exist
        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            self.assertEqual(d.get_sci_prepped_lig(temp_dir,d.SCI_PREPPED_LIG_SUFFIX), False)
        finally:
            shutil.rmtree(temp_dir)    

        # Test where correct file exists
        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            ligfile = os.path.join(temp_dir, 'lig_TST'+d.SCI_PREPPED_LIG_SUFFIX)
            f = open(ligfile, 'w')
            f.write('There is no test to ensure the contents are valid. This string should pass though.')
            f.flush()
            f.close()
            self.assertEqual(d.get_sci_prepped_lig(temp_dir,d.SCI_PREPPED_LIG_SUFFIX), ligfile)
        finally:
            shutil.rmtree(temp_dir)    


        # Test where file has wrong suffix
        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            ligfile = os.path.join(temp_dir, 'lig_TST_prepared.xxx')
            f = open(ligfile, 'w')
            f.write('There is no test to ensure the contents are valid. This string should pass though.')
            f.flush()
            f.close()
            self.assertEqual(d.get_sci_prepped_lig(temp_dir,d.SCI_PREPPED_LIG_SUFFIX), False)
        finally:
            shutil.rmtree(temp_dir)    


        # Test where multiple correct files exist
        temp_dir = tempfile.mkdtemp()
        try:
            d = Dock()
            ligfile = os.path.join(temp_dir, 'lig_TST'+d.SCI_PREPPED_LIG_SUFFIX)
            f = open(ligfile, 'w')
            f.write('There is no test to ensure the contents are valid. This string should pass though.')
            f.flush()
            f.close()

            ligfile = os.path.join(temp_dir, 'lig_TS2'+d.SCI_PREPPED_LIG_SUFFIX)
            f = open(ligfile, 'w')
            f.write('There is no test to ensure the contents are valid. This string should pass though.')
            f.flush()
            f.close()
            self.assertEqual(d.get_sci_prepped_lig(temp_dir,d.SCI_PREPPED_LIG_SUFFIX), False)
        finally:
            shutil.rmtree(temp_dir)    


    def test_parse_lig_filename(self):
        d = Dock()

        # Blank name
        self.assertEqual(d.parse_lig_filename(''), False)

        # Good name
        ligfile = os.path.join('some_directory', 'lig_TST'+d.SCI_PREPPED_LIG_SUFFIX)
        self.assertEqual(d.parse_lig_filename(ligfile), 'TST')

        # Good name, but in list
        self.assertEqual(d.parse_lig_filename([ligfile]), False)

        # Wrong suffix
        ligfile = os.path.join('some_directory', 'lig_TST_prepared.xxx')
        self.assertEqual(d.parse_lig_filename(ligfile), False)

        # Change expected suffix to "_prepared.xxx"
        d.SCI_PREPPED_LIG_SUFFIX = '_prepared.xxx'
        self.assertEqual(d.parse_lig_filename(ligfile), 'TST')


    def test_parse_cand_name(self):
        d = Dock()
        
        # Blank name
        self.assertEqual(d.parse_cand_name(''), False)
        
        # Good name
        candfile = os.path.join('some_directory', 'hiResHolo-1fcz_1fcy'+d.SCI_PREPPED_PROT_SUFFIX)
        self.assertEqual(d.parse_cand_name(candfile), ('hiResHolo','1fcz','1fcy'))
        
        # Good name, but in list
        self.assertEqual(d.parse_cand_name([candfile]), False)
        
        # Wrong suffix
        candfile = os.path.join('some_directory', 'hiResHolo-1fcz_1fcy_prepared.xxx')
        self.assertEqual(d.parse_cand_name(candfile), False)
        
        # Change to be right suffix
        d.SCI_PREPPED_PROT_SUFFIX = '_prepared.xxx'
        self.assertEqual(d.parse_cand_name(candfile), ('hiResHolo','1fcz','1fcy'))

        
    def test_run_dock(self):
        orig_dir = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        temp_dir_abs = os.path.abspath(temp_dir)
        prot_sci_prep_dir = os.path.join(temp_dir,'prot_sci_prep_dir')
        lig_sci_prep_dir = os.path.join(temp_dir,'lig_sci_prep_dir')
        dock_dir = os.path.join(temp_dir,'dock_dir')
        try:        
            os.mkdir(prot_sci_prep_dir)
            sample_prot_dir = os.path.join(prot_sci_prep_dir,'1fcz')
            os.mkdir(sample_prot_dir)
            sample_prot_file = os.path.join(sample_prot_dir, 'hiResHolo-1fcz_1fcy_prepared.pdb')
            with open(sample_prot_file, 'wb') as of:
                of.write('test')
            sample_center_file = os.path.join(sample_prot_dir, 'center.txt')
            with open(sample_center_file, 'wb') as of:
                of.write('1,2,3')
            sample_cand_targinfo_file = os.path.join(sample_prot_dir, '1fcz.txt')
            with open(sample_cand_targinfo_file, 'wb') as of:
                of.write('query, 1fcz')
        
            os.mkdir(lig_sci_prep_dir)
            sample_lig_dir = os.path.join(lig_sci_prep_dir,'1fcz')
            os.mkdir(sample_lig_dir)
            sample_lig_file = os.path.join(sample_lig_dir, 'lig_156_prepared.sdf')
            with open(sample_lig_file, 'wb') as of:
                of.write('test')
            sample_lig_targinfo_file = os.path.join(sample_lig_dir, '1fcz.txt')
            with open(sample_lig_targinfo_file, 'wb') as of:
                of.write('query, 1fcz')            
            os.mkdir(dock_dir)
            d = Dock()
            d.run_dock(prot_sci_prep_dir, lig_sci_prep_dir, dock_dir)
            
            dock_targ_dir = os.path.join(dock_dir,'1fcz')
            self.assertEqual(os.path.isdir(dock_targ_dir), True)
            dock_cand_prep_dir = os.path.join(dock_targ_dir,'hiResHolo_1fcy_tech_prep')
            self.assertEqual(os.path.isdir(dock_cand_prep_dir), True)
            dock_lig_prep_dir = os.path.join(dock_targ_dir,'lig_156_tech_prep')
            self.assertEqual(os.path.isdir(dock_lig_prep_dir), True)
            dock_cand_dock_dir = os.path.join(dock_targ_dir,'hiResHolo_1fcy_docking')
            self.assertEqual(os.path.isdir(dock_cand_dock_dir), True)
            dock_lig_file = os.path.join(dock_cand_dock_dir,'lig_156_prepared.sdf')
            self.assertEqual(os.path.exists(dock_lig_file), True)
            dock_cand_file = os.path.join(dock_cand_dock_dir,'hiResHolo-1fcz_1fcy_prepared.pdb')
            self.assertEqual(os.path.exists(dock_cand_file), True)
        finally:
            os.chdir(orig_dir)
            shutil.rmtree(temp_dir_abs)

   
           
    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
