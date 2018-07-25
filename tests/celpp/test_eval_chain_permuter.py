#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'jwagner'


import unittest
import tempfile
import os.path
import shutil
import os
import sys

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.blastnfilter import BlastNFilterTask
from d3r.celpp.challengedata import ChallengeDataTask
from d3r.celpp.evaluation import EvaluationTaskFactory
from d3r.celpp.evaluation import EvaluationTask
from d3r.celpp.participant import ParticipantDatabase
from d3r.celpp.participant import Participant
from d3r.celpp.task import SmtpEmailer
from d3r.celpp.evaluation import EvaluationEmailer
from mock import Mock


class TestChainPermuter(unittest.TestCase):
    def setUp(self):
        try:
            import openeye.oechem  # noqa: F401
            sys.stdout.write("openeye.oechem importable.\n")
            self.has_oechem = True
        except Exception:
            sys.stdout.write("openeye.oechem not importable.\n")
            self.has_oechem = False

        ret_code = os.system('PYTHONPATH= $SCHRODINGER/run '
                             'split_structure.py &> '
                             '/dev/null')
        if ret_code == 512:
            sys.stdout.write("Schrodinger commandline available.\n")
            self.has_schrodinger = True
        else:
            sys.stdout.write("Schrodinger commandline not available.\n")
            self.has_schrodinger = False

    def test_permute_chains_and_take_min_rmsd(self):
        if not(self.has_oechem) or not (self.has_schrodinger):
            sys.stdout.write("Missing essential component. Skipping test.\n")
            return
        test_scale = 'github'
        method_list = ['autodockvina']
        try:
            temp_dir = tempfile.mkdtemp()
            for method in method_list:
                params = D3RParameters()
                params.evaluation = 'evaluate.py'
                params.pdbdb = os.path.abspath('tests/celpp/eval_test_data/'
                                               '%s_test_data/mini_pdb/'
                                               % test_scale)

                # Make blastnfilter task
                blastnfiltertask = BlastNFilterTask(temp_dir, params)
                blastnfiltertask.create_dir()
                open(os.path.join(blastnfiltertask.get_dir(),
                                  D3RTask.COMPLETE_FILE), 'a').close()
                source_dir = os.path.abspath('tests/celpp/eval_test_data/'
                                             '%s_test_data/stage.3.'
                                             'blastnfilter/'
                                             % test_scale)
                os.system('cp -r %s/* %s' % (source_dir,
                                             blastnfiltertask.get_dir()))

                # Make challengedata task
                challengedatatask = ChallengeDataTask(temp_dir, params)
                challengedatatask.create_dir()
                open(os.path.join(challengedatatask.get_dir(),
                                  D3RTask.COMPLETE_FILE), 'a').close()
                source_dir = os.path.abspath('tests/celpp/eval_test_data/'
                                             '%s_test_data/stage.4.'
                                             'challengedata/'
                                             % test_scale)
                os.system('cp -r %s/* %s' % (source_dir,
                                             challengedatatask.get_dir()))

                # Make dock task
                docktask = D3RTask(temp_dir, params)
                docktask.set_name(method)
                docktask.set_stage(EvaluationTaskFactory.DOCKSTAGE)
                docktask.create_dir()
                open(os.path.join(docktask.get_dir(), D3RTask.COMPLETE_FILE),
                     'a').close()

                source_dir = os.path.abspath('tests/celpp/eval_test_data/'
                                             '%s_test_data/stage.6.%s/'
                                             % (test_scale, method))
                os.system('cp -r %s/* %s' % (source_dir, docktask.get_dir()))
                print 'docktask.get_dir()', docktask.get_dir()
                print 'temp_dir', temp_dir
                task = EvaluationTask(temp_dir, '%s.evaluation' % (method),
                                      docktask, params)

                plist = [Participant('1name', '1d3rusername', '12345',
                                     'bob@bob.com,joe@joe.com')]
                smtpemailer = SmtpEmailer()
                mockserver = D3RParameters()
                mockserver.sendmail = Mock()
                mockserver.quit = Mock()
                smtpemailer.set_alternate_smtp_server(mockserver)
                emailer = EvaluationEmailer(ParticipantDatabase(plist), None)
                emailer.set_alternate_smtp_emailer(smtpemailer)
                task.set_evaluation_emailer(emailer)

                task.run()
                val = task.get_evaluation_summary()
                self.assertEqual(val,
                                 '\nEvaluation of docking\n==================='
                                 '==\nTarge'
                                 't_PDBID        LMCSS            SMCSS       '
                                 '     h'
                                 'iResApo         hiResHolo        hiTanimoto '
                                 '      '
                                 'LMCSS_ori_distance \n\nSummary Statistics\n'
                                 '\nNumber_of'
                                 '_cases     0                0               '
                                 ' 1    '
                                 '            0                0              '
                                 '  \nAve'
                                 'rage                                        '
                                 '      '
                                 ' 6.447                                     '
                                 '       '
                                 '  \nMaximum                                '
                                 '        '
                                 '       6.447                               '
                                 '       '
                                 '        \nMinimum                           '
                                 '       '
                                 '             6.447                          '
                                 '      '
                                 '              \nMedian                      '
                                 '       '
                                 '                   6.447                    '
                                 '      '
                                 '                    \n\nIndividual Results\n'
                                 '\n5t6d    '
                                 '                                            '
                                 '  6.44'
                                 '7 (2.176 )                                  '
                                 '      '
                                 '              \n\n')
                self.assertEqual(task.get_error(), None)
        finally:
            # pass
            shutil.rmtree(temp_dir)


if __name__ == '__main__':
    unittest.main()
