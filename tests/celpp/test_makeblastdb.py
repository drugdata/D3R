#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'

import os.path
import tempfile
import unittest
import shutil
import gzip
import itertools
import string

"""
test_makeblastdb
--------------------------------

Tests for `makeblastdb` module.
"""

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import D3RTask
from d3r.celpp.makeblastdb import MakeBlastDBTask
from tests.celpp import test_task


class TestMakeBlastDBTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_constructor(self):
        params = D3RParameters()
        task = MakeBlastDBTask('/foo', params)
        self.assertEqual(task.get_name(), 'makeblastdb')
        self.assertEqual(task.get_stage(), 1)
        self.assertEqual(task.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(task.get_path(), '/foo')
        self.assertEqual(task.get_dir_name(), 'stage.1.makeblastdb')
        test_task.try_update_status_from_filesystem(self, task)

    def test_get_uploadable_files(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()

            # test on empty dir
            self.assertEqual(task.get_uploadable_files(), [])

            # test with stderr/stdout files
            stdout = os.path.join(task.get_dir(), 'makeblastdb.stdout')
            open(stdout, 'a').close()
            stderr = os.path.join(task.get_dir(), 'makeblastdb.stderr')
            open(stderr, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(stdout)
            flist.index(stderr)

            # test with pdb_seqres.txt.gz file
            open(task.get_pdb_seqres_txt_gz(), 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(stdout)
            flist.index(stderr)
            flist.index(task.get_pdb_seqres_txt_gz())

            os.unlink(task.get_pdb_seqres_txt_gz())

            # test where pdb_seqres.txt is a dir unlikely but why not check
            os.makedirs(task.get_pdb_seqres_txt_gz())
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 2)
            flist.index(stdout)
            flist.index(stderr)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_pdbseqres_txt_gz(self):
        params = D3RParameters()
        task = MakeBlastDBTask('/foo', params)
        self.assertEqual(task.get_pdb_seqres_txt_gz(),
                         os.path.join('/foo', task.get_dir_name(),
                                      'pdb_seqres.txt.gz'))

    def test_get_pdbseqres_txt(self):
        params = D3RParameters()
        task = MakeBlastDBTask('/foo', params)
        self.assertEqual(task.get_pdb_seqres_txt(),
                         os.path.join('/foo', task.get_dir_name(),
                                      'pdb_seqres.txt'))

    def test_get_sequence_count_no_file(self):
        params = D3RParameters()
        task = MakeBlastDBTask('/foo', params)

        self.assertEqual(task._get_sequence_count_message(),
                         '# sequence(s): Error unable to parse file')

    def test_get_sequence_count_file_has_zero_size(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            open(task.get_pdb_seqres_txt(), 'a').close()
            self.assertEqual(task._get_sequence_count_message(),
                             '# sequence(s): 0')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_sequence_count_file_has_no_seqs(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            f.write('hi\nhow\nare\n')
            f.flush()
            f.close()
            self.assertEqual(task._get_sequence_count_message(),
                             '# sequence(s): 0')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_sequence_count_file_has_one_seq(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            f.write('hi\n>aseq\nare\n')
            f.flush()
            f.close()
            self.assertEqual(task._get_sequence_count_message(),
                             '# sequence(s): 1')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_sequence_count_file_has_multiple_seqs(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            f.write('>hi\n>seq\n>are\n')
            f.flush()
            f.close()
            self.assertEqual(task._get_sequence_count_message(),
                             '# sequence(s): 3')
        finally:
            shutil.rmtree(temp_dir)

    def test_get_sequence_count_file_no_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            pdbid = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbid), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pbdid_from_pdb_seqres_txt_no_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            pdbset = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pbdid_from_pdb_seqres_txt_empty_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            open(task.get_pdb_seqres_txt(), 'a').close()
            self.assertEqual(os.path.isfile(task.get_pdb_seqres_txt()), True)
            pdbset = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pbdid_from_pdb_seqres_txt_file_no_seqs(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            f.write('hi\nhow\nare\nyou')
            f.flush()
            f.close()
            self.assertEqual(os.path.isfile(task.get_pdb_seqres_txt()), True)
            pdbset = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pbdid_from_pdb_seqres_txt_with_seqs(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            f.write('>101m_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVK'
                    'HLKTEAEMKASEDLKKHG\n')
            f.write('>102l_A mol:protein length:165  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKA'
                    'IGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAL'
                    'INMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRV'
                    'ITTFRTGTWDAYKNL\n')
            f.write('>102l_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHL'
                    'KTEAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKI'
                    'PIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG\n')
            f.write('>103l_A mol:protein length:167  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNSLDAAKSELD'
                    'KAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRA'
                    'ALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAK'
                    'RVITTFRTGTWDAYKNL\n')
            f.write('>10jj3m_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRF\n')
            f.write('>104l_A mol:protein length:166  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSAAE\n')
            f.write('>104l_B mol:protein length:166  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSAKNL\n')

            f.flush()
            f.close()
            self.assertEqual(os.path.isfile(task.get_pdb_seqres_txt()), True)
            pdbset = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbset), 4)
            self.assertEqual('101M' in pdbset, True)
            self.assertEqual('102L' in pdbset, True)
            self.assertEqual('103L' in pdbset, True)
            self.assertEqual('104L' in pdbset, True)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pbdid_from_pdb_seqres_txt_wrong_len_pdbids(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            f.write('>1m_A mol:protein length:154  MYOGLOBIN\n')
            f.write('MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVK'
                    'HLKTEAEMKASEDLKKHG\n')
            f.write('>abcdel_A mol:protein length:165  T4 LYSOZYME\n')
            f.write('MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKA'
                    'IGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAL'
                    'INMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRV'
                    'ITTFRTGTWDAYKNL\n')

            f.flush()
            f.close()
            self.assertEqual(os.path.isfile(task.get_pdb_seqres_txt()), True)
            pdbset = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbset), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_get_set_of_pbdid_from_pdb_seqres_txt_with_400k_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            f = open(task.get_pdb_seqres_txt(), 'w')
            perms = itertools.permutations(string.ascii_lowercase +
                                           '123456789', 4)
            counter = 0
            try:
                limit = 400000
                while counter < limit:
                    f.write('>' + ''.join(map(str, perms.next())) +
                            '_A mol:protein length:165  T4 LYSOZYME\n')
                    f.write('MVLSEGEWQLVLH\n')
                    counter += 1
            except StopIteration:
                pass
            f.flush()
            f.close()

            pdbset = task.get_set_of_pbdid_from_pdb_seqres_txt()
            self.assertEqual(len(pdbset), counter)

        finally:
            shutil.rmtree(temp_dir)

    def test_can_run_where_task_is_complete(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            open(os.path.join(task.get_dir(), 'complete'), 'a').close()
            self.assertEqual(task.can_run(), False)
        finally:
            shutil.rmtree(temp_dir)

    def test_can_run_where_task_failed(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            open(os.path.join(task.get_dir(), 'error'), 'a').close()
            self.assertEqual(task.can_run(), False)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_can_run_is_false(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.create_dir()
            open(os.path.join(task.get_dir(), 'error'), 'a').close()
            task.run()

            self.assertEqual(task._can_run, False)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_pdbsequrl_is_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = MakeBlastDBTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(), 'cannot download files cause '
                                               'pdbsequrl not set')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_makeblastdb_is_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbsequrl = 'pdbsequrl'
            task = MakeBlastDBTask(temp_dir, params)
            task.run()
            self.assertEqual(task.get_error(), 'cannot make blast database '
                                               'cause makeblastdb not set')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_download_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.pdbsequrl = 'file://doesnotexist'
            params.makeblastdb = 'makeblastdb'
            task = MakeBlastDBTask(temp_dir, params)
            task._retrysleep = 0
            task._maxretries = 1
            task.run()
            self.assertEqual(task.get_error(), 'Unable to download file: ' +
                             'file://doesnotexist')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_gunzip_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            fakegz = os.path.join(temp_dir, 'fake.gz')

            f = open(fakegz, 'w')
            f.write('hello\n')
            f.flush()
            f.close()

            params.pdbsequrl = 'file://'+fakegz
            params.makeblastdb = 'makeblastdb'
            task = MakeBlastDBTask(temp_dir, params)
            task._retrysleep = 0
            task._maxretries = 1
            task.run()
            self.assertEqual(task.get_error(), 'Unable to uncompress file: ' +
                             task.get_pdb_seqres_txt())
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_makeblastdb_fails(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            fakegz = os.path.join(temp_dir, 'fake.gz')

            f = gzip.open(fakegz, 'wb')
            f.write('hello\n')
            f.flush()
            f.close()

            params.pdbsequrl = 'file://'+fakegz
            params.makeblastdb = 'false'
            task = MakeBlastDBTask(temp_dir, params)
            task._retrysleep = 0
            task._maxretries = 1
            task.run()
            self.assertEqual(task.get_error(), 'Non zero exit code: 1 '
                                               'received. Standard out:'
                                               '  Standard error: ')
        finally:
            shutil.rmtree(temp_dir)

    def test_run_where_everything_is_successful(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            fakegz = os.path.join(temp_dir, 'fake.gz')

            f = gzip.open(fakegz, 'wb')
            f.write('hello\n')
            f.flush()
            f.close()

            params.pdbsequrl = 'file://'+fakegz
            params.makeblastdb = 'echo'
            task = MakeBlastDBTask(temp_dir, params)
            task._retrysleep = 0
            task._maxretries = 1
            task.run()
            self.assertEqual(task.get_error(), None)

            # check echo.stdout file for valid arguments
            f = open(os.path.join(task.get_dir(), 'echo.stdout'), 'r')
            line = f.readline()

            self.assertEqual(line, '-in ' +
                             task.get_pdb_seqres_txt() +
                             ' -out ' +
                             os.path.join(task.get_dir(), 'pdb_db') +
                             ' -dbtype prot\n')

            f.close()

            lines = task.get_email_log().split('\n')
            self.assertEqual(lines[2], '# sequence(s): 0')
            f.close()
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
