__author__ = 'j5wagner'

import unittest
import tempfile
import os.path

"""
test_evaluation
--------------------------------

Tests for `evaluation` module.
"""
import shutil
import os

#from mock import Mock

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
#from d3r.celpp.dataimport import DataImportTask
#from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.evaluation import EvaluationTask
#from d3r.celpp.glide import GlideTask
#from d3r.celpp.evaluation import PathNotDirectoryError
#from d3r.celpp.participant import ParticipantDatabaseFromCSVFactory
#from d3r.celpp.participant import ParticipantDatabase
#from d3r.celpp.participant import Participant
#from d3r.celpp.task import SmtpEmailer
#from d3r.celpp.evaluation import EvaluationEmailer


class TestChainPermuter(unittest.TestCase):
   def setUp(self):
      try:
         import openeye.oechem
         print "epeneye.oechem importable."
         self.has_oechem = True
      except:
         print "epeneye.oechem not importable."
         self.has_oechem = False
      ret_code = os.system('$SCHRODINGER/run split_structure.py &> /dev/null')
      if ret_code == 512:
         print "Schrodinger commandline available."
         self.has_schrodinger = True
      else:
         print "Schrodinger commandline not available."
         self.has_schrodinger = False



   #@unittest.skipIf(not(self.hasOechem), "openeye.oechem not available. "
   #                                      "Skipping eval performance test")
   def test_permute_chains_and_take_min_rmsd(self):
      if not(self.has_oechem) or not (self.has_schrodinger):
         print "Missing essential component. Skipping test."
         return
      #test_scale = 'minimal'
      #method_list = ['33567.extsubmission']
      test_scale = 'github'
      method_list = ['autodockvina']
      #test_scale = 'maximal'
      #method_list = ['33674.extsubmission','33674_2.extsubmission']
      #method_list = ['33567.extsubmission','33567_2.extsubmission','33674.extsubmission','33674_2.extsubmission','glide','autodockvina']
      try:
         temp_dir = tempfile.mkdtemp()
         for method in method_list:
            #print "Evaluating", method
            params = D3RParameters()
            params.evaluation = 'evaluate.py'
            params.pdbdb = os.path.abspath('tests/celpp/eval_test_data/%s_test_data/mini_pdb/' %(test_scale))
            docktask = D3RTask(temp_dir, params)
            docktask.set_name(method)
            docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
            docktask.create_dir()
            open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                 'a').close()

            source_dir = os.path.abspath('tests/celpp/eval_test_data/%s_test_data/stage.6.%s/' %(test_scale, method))
            os.system('cp -r %s/* %s' %(source_dir, docktask.get_dir()))
            task = EvaluationTask(temp_dir, 
                                  '%s.evaluation' %(method), 
                                  docktask, 
                                  params)
            #print dir(task)
            task.run()
            val = task.get_evaluation_summary()
            self.assertEqual(val, '\nEvaluation of docking\n=====================\nTarget_PDBID        LMCSS     SMCSS     hiResApo  hiResHolo hiTanimoto \nNumber_of_cases     0         0         1         0         0         \nAverage                                 6.452                         \nMinimum                                 6.452                         \nMaximum                                 6.452                         \n5t6d                                    6.452                         \n\n')
            self.assertEqual(task.get_error(),None)
      finally:
         shutil.rmtree(temp_dir)
         #print temp_dir
         print 'done'

if __name__ == '__main__':
    unittest.main()
