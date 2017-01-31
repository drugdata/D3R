#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'churas'


import unittest
import os
import tempfile
import shutil


"""
test_participant
----------------------------------

Tests for `participant` module.
"""

from d3r.celpp.participant import Participant
from d3r.celpp.participant import ParticipantDatabase
from d3r.celpp.participant import ParticipantDatabaseFromCSVFactory


class TestParticipant(unittest.TestCase):
    def setUp(self):
        pass

    def test_participant(self):
        p = Participant(None, None, None, None)
        self.assertEqual(p.get_d3rusername(), None)
        self.assertEqual(p.get_email(), None)
        self.assertEqual(p.get_guid(), None)
        self.assertEqual(p.get_name(), None)

        p = Participant('name', 'd3rusername', 'guid', 'email@email.com')
        self.assertEqual(p.get_d3rusername(), 'd3rusername')
        self.assertEqual(p.get_email(), 'email@email.com')
        self.assertEqual(p.get_guid(), 'guid')
        self.assertEqual(p.get_name(), 'name')

    def test_participant_database(self):

        pdb = ParticipantDatabase(None)
        self.assertEqual(pdb.get_participants(), None)
        self.assertEqual(pdb.get_participant_by_guid('bob'), None)

        pdb = ParticipantDatabase([])
        self.assertEqual(len(pdb.get_participants()), 0)
        self.assertEqual(pdb.get_participant_by_guid('bob'), None)
        self.assertEqual(pdb.get_participant_by_guid('bob',
                                                     exact_match=True), None)
        self.assertEqual(pdb.get_participant_by_guid('guidx'), None)
        self.assertEqual(pdb.get_participant_by_guid('xguid'), None)
        self.assertEqual(pdb.get_participant_by_guid('GUID'), None)

        plist = [Participant('1name', '1d3rusername', 'guid',
                             '1email@email.com'),
                 Participant('2name', '2d3rusername', 'guidy',
                             '2email@email.com'),
                 Participant('3name', '3d3rusername', '12345',
                             '3email@email.com')]

        pdb = ParticipantDatabase(plist)
        self.assertEqual(pdb.get_participant_by_guid(None), None)
        self.assertEqual(pdb.get_participant_by_guid('bob'), None)
        self.assertEqual(pdb.get_participant_by_guid('guidx'), None)
        self.assertEqual(pdb.get_participant_by_guid('xguid'), None)
        self.assertEqual(pdb.get_participant_by_guid('GUID'), None)

        p = pdb.get_participant_by_guid('guid')
        self.assertEqual(p.get_name(), '1name')

        p = pdb.get_participant_by_guid('guidy')
        self.assertEqual(p.get_name(), '2name')

        p = pdb.get_participant_by_guid('12345')
        self.assertEqual(p.get_name(), '3name')

        p = pdb.get_participant_by_guid('12345', exact_match=True)
        self.assertEqual(p.get_name(), '3name')

        p = pdb.get_participant_by_guid('12345_2')
        self.assertEqual(p.get_name(), '3name')

        p = pdb.get_participant_by_guid('12345_2', exact_match=False)
        self.assertEqual(p.get_name(), '3name')

        p = pdb.get_participant_by_guid('12345_mydock', exact_match=False)
        self.assertEqual(p.get_name(), '3name')

        p = pdb.get_participant_by_guid('12345_my-dock', exact_match=False)
        self.assertEqual(p.get_name(), '3name')

        # Suffixes can't contain underscores
        p = pdb.get_participant_by_guid('12345_my_dock', exact_match=False)
        self.assertEqual(p, None)

        p = pdb.get_participant_by_guid('12345_2', exact_match=True)
        self.assertEqual(p, None)

        p = pdb.get_participant_by_guid('12345_2_b')
        self.assertEqual(p, None)

    def test_participant_database_from_csv_factory(self):
        temp_dir = tempfile.mkdtemp()
        try:

            csvfile = os.path.join(temp_dir, 'participant_list.csv')

            # parse None as file
            pfac = ParticipantDatabaseFromCSVFactory(None)
            self.assertEqual(pfac.get_participant_database(), None)

            # parse non existant file
            pfac = ParticipantDatabaseFromCSVFactory(csvfile)
            self.assertEqual(pfac.get_participant_database(), None)

            # parse empty file
            open(csvfile, 'a').close()
            pfac = ParticipantDatabaseFromCSVFactory(csvfile)
            pdb = pfac.get_participant_database()
            self.assertEqual(len(pdb.get_participants()), 0)

            # parse file with only header
            f = open(csvfile, 'w')
            f.write('name,d3rusername,guid,email\n')
            f.flush()
            f.close()
            pfac = ParticipantDatabaseFromCSVFactory(csvfile)
            pdb = pfac.get_participant_database()
            self.assertEqual(len(pdb.get_participants()), 0)

            # parse file with header and 1 valid entry
            f = open(csvfile, 'w')
            f.write('name,d3rusername,guid,email\n')
            f.write('joe,jj,123,j@j.com\n')
            f.flush()
            f.close()
            pfac = ParticipantDatabaseFromCSVFactory(csvfile)
            pdb = pfac.get_participant_database()
            self.assertEqual(len(pdb.get_participants()), 1)
            p = pdb.get_participant_by_guid('123')
            self.assertEqual(p.get_email(), 'j@j.com')

            # parse file with no header and 1 valid entry
            f = open(csvfile, 'w')
            f.write('joe,jj,123,j@j.com\n')
            f.flush()
            f.close()
            pfac = ParticipantDatabaseFromCSVFactory(csvfile)
            pdb = pfac.get_participant_database()
            self.assertEqual(len(pdb.get_participants()), 1)
            p = pdb.get_participant_by_guid('123')
            self.assertEqual(p.get_email(), 'j@j.com')

            # parse file with header and 2 valid and 1 invalid entries

            f = open(csvfile, 'w')
            f.write('name,d3rusername,guid,email\n')
            f.write('joe,jj,123,j@j.com\n')
            f.write('uhoh,2,many,commas,uhoh@email.com\n')
            f.write('phil,pp,456,p@p.com\n')

            f.flush()
            f.close()
            pfac = ParticipantDatabaseFromCSVFactory(csvfile)
            pdb = pfac.get_participant_database()
            self.assertEqual(len(pdb.get_participants()), 2)
            p = pdb.get_participant_by_guid('123')
            self.assertEqual(p.get_email(), 'j@j.com')
            p = pdb.get_participant_by_guid('456')
            self.assertEqual(p.get_email(), 'p@p.com')

            # parse file with carriage returns
            f = open(csvfile, 'w')
            f.write('name,d3rusername,guid,email\r')
            f.write('joe bob , jb,12345,j@j.com\r')
            f.write('J W,ha@ha.com,33567 , ha@ha.com')
            f.flush()
            f.close()
            pdb = pfac.get_participant_database()
            self.assertEqual(len(pdb.get_participants()), 2)
            p = pdb.get_participant_by_guid('12345')
            self.assertEqual(p.get_email(), 'j@j.com')
            self.assertEqual(p.get_name(), 'joe bob')
            self.assertEqual(p.get_d3rusername(), 'jb')
            p = pdb.get_participant_by_guid('33567')
            self.assertEqual(p.get_email(), 'ha@ha.com')
            self.assertEqual(p.get_name(), 'J W')
            self.assertEqual(p.get_d3rusername(), 'ha@ha.com')

        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
