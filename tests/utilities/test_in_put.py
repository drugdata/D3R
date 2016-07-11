__author__ = 'churas'

import unittest
import tempfile
import os.path

from d3r.blast.query import Query
"""
test_in_put
--------------------------------

Tests for `in_put` module.
"""

import shutil
from d3r.utilities import in_put


class TestInPut(unittest.TestCase):
    def setUp(self):
        pass

    def get_valid_new_release_sequence_tsv(self):
        val = """PDB_ID  Sequence_Count  Sequence
2N1I    1       GSGFPTSEDFTPKEGSPYEAPVYIPEDIPIPADFELRESSIPGAGLGVWAKRKMEAG
2N27    1       ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADEMIR
2N3M    1       (DT)(DG)(DG)(DT)(DG)(DG)(DT)(DG)(DG)(DT)(DT)(DG)(DT)(DT)(DG)
2N4M    1       (DG)(DT)(DG)(DC)(4E9)(DT)(DG)(DT)(DT)(DT)(DG)(DT)
2N4M    2       (DA)(DC)(DA)(DA)(DA)(DC)(DA)(DC)(DG)(DC)(DA)(DC)
2N7A    1       MEGDRQYGDGYLLQVQELVTVQEGLSVHVPCSFSYPQDGWTDSDPVHGYWFRAGDRPYQDAP
2N7B    1       MEGDRQYGDGYLLQVQELVTVQEGLSVHVPCSFSYPQDGWTDSDPVHGYWFRAGDRPYQDAP
2N9T    1       YCQKWMWTCDSERKCCEGMVCRLWCKKKLW"""
        return val

    def test_read_sequences_valid_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tsvdata = self.get_valid_new_release_sequence_tsv()
            poly = os.path.join(temp_dir, 'foo.tsv')
            f = open(poly, 'w')
            f.write(tsvdata)
            f.flush()
            f.close()
            queries = in_put.read_sequences(poly)
            self.assertEqual(len(queries), 7)
            self.assertEqual(queries[0].pdb_id, '2n1i')
            self.assertEqual(queries[6].pdb_id, '2n9t')
        finally:
            shutil.rmtree(temp_dir)

    def test_read_sequences_invalid_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            poly = os.path.join(temp_dir, 'foo.tsv')
            f = open(poly, 'w')
            f.write('blah\nblah\n')
            f.flush()
            f.close()
            queries = in_put.read_sequences(poly)
            self.assertEqual(len(queries), 0)
        finally:
            shutil.rmtree(temp_dir)

    def get_valid_nonpolymer_tsv(self):
        val = """PDB_ID  Component_ID    InChI
2N1I    CA      InChI=1S/Ca/q+2
2N27    4DY     InChI=1S/C18H27NO3/c1-14(2)8-6-4-5-7-9-18(21)19-13-15-10
2N7B    0D8     InChI=1S/C3H9NO/c4-2-1-3-5/h5H,1-4H2
4XGS    OFO     InChI=1S/2Fe.H2O.O/h;;1H2;/q;+1;;/p-1
4YDE    EDO     InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2
4YDE    ZN      InChI=1S/Zn/q+2"""
        return val

    def test_read_ligands_valid_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tsvdata = self.get_valid_nonpolymer_tsv()
            nonpoly = os.path.join(temp_dir, 'foo2.tsv')
            f = open(nonpoly, 'w')
            f.write(tsvdata)
            f.flush()
            f.close()
            queries = []
            q = Query()
            q.pdb_id = '2n1i'
            queries.append(q)
            q = Query()
            q.pdb_id = '2n27'
            queries.append(q)
            q = Query()
            q.pdb_id = '2n7b'
            queries.append(q)
            in_put.read_ligands(nonpoly, queries)

        finally:
            shutil.rmtree(temp_dir)

    def get_valid_ph(self):
        val = """PDB_ID  _exptl_crystal_grow.pH
2N1I    6.5
2N27    8
2N7B    8
2N3M    5.5
2N4M    8.5
2N7A    8"""
        return val

    def test_read_ph_valid_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tsvdata = self.get_valid_ph()
            ph = os.path.join(temp_dir, 'foo.tsv')
            f = open(ph, 'w')
            f.write(tsvdata)
            f.flush()
            f.close()
            queries = []
            q = Query()
            q.pdb_id = '2n1i'
            queries.append(q)
            q = Query()
            q.pdb_id = '2n27'
            queries.append(q)
            q = Query()
            q.pdb_id = '2n7b'
            queries.append(q)

            in_put.read_ph(ph, queries)

            self.assertEqual(queries[0].exp_ph, '6.5')
            self.assertEqual(queries[1].exp_ph, '8')
            self.assertEqual(queries[2].exp_ph, '8')

        finally:
            shutil.rmtree(temp_dir)

    def test_create_queries(self):
        temp_dir = tempfile.mkdtemp()
        try:
            poly = self.get_valid_new_release_sequence_tsv()
            poly_file = os.path.join(temp_dir, 'poly.tsv')
            f = open(poly_file, 'w')
            f.write(poly)
            f.flush()
            f.close()

            nopoly = self.get_valid_nonpolymer_tsv()
            nopoly_file = os.path.join(temp_dir, 'nopoly.tsv')
            f = open(nopoly_file, 'w')
            f.write(nopoly)
            f.flush()
            f.close()

            ph = self.get_valid_ph()
            ph_file = os.path.join(temp_dir, 'ph.tsv')
            f = open(ph_file, 'w')
            f.write(ph)
            f.flush()
            f.close()

            queries = in_put.create_queries(poly_file,
                                            nopoly_file, ph_file)
            self.assertEqual(len(queries), 7)
            self.assertEqual(queries[0].pdb_id, '2n1i')
            self.assertEqual(queries[6].pdb_id, '2n9t')

        finally:
            shutil.rmtree(temp_dir)

    def test_label(self):
        self.assertEqual(in_put.label('VA3'), 'do_not_call')
        self.assertEqual(in_put.label('hahaha'), 'dock')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
