#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import logging
import os.path

import d3r.task

"""
test_task
----------------------------------

Tests for `task` module.
"""

import unittest

from d3r import task


class TestD3rTask(unittest.TestCase):

    def setUp(self):
        pass

    def test_find_latest_year(self):
        tempDir = tempfile.mkdtemp()

	
        self.assertEqual(task._find_latest_year(tempDir),None)
       
        os.mkdir(os.path.join(tempDir,"foo"))
        os.mkdir(os.path.join(tempDir,"2014"))
        self.assertEqual(task._find_latest_year(tempDir),
                         os.path.join(tempDir,'2014')) 

        open(os.path.join(tempDir,'2016'),'a').close()
        self.assertEqual(task._find_latest_year(tempDir),  
                         os.path.join(tempDir,'2014'))

        os.mkdir(os.path.join(tempDir,'2012'))
           
        self.assertEqual(task._find_latest_year(tempDir),
                         os.path.join(tempDir,'2014'))

        os.mkdir(os.path.join(tempDir,'2015'))


        self.assertEqual(task._find_latest_year(tempDir),
                         os.path.join(tempDir,'2015'))

        os.rmdir(os.path.join(tempDir,"foo"))
        os.rmdir(os.path.join(tempDir,"2012"))
        os.rmdir(os.path.join(tempDir,"2014"))
        os.rmdir(os.path.join(tempDir,"2015"))
        os.remove(os.path.join(tempDir,'2016'))
        
        os.rmdir(tempDir)
     
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
