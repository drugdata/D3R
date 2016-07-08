__author__ = 'churas'

import unittest
import tempfile
import os.path

"""
test_blastnfilter
--------------------------------

Tests for `blastnfilter` module.
"""

import shutil
from d3r.utilities import in_put


class TestBlastNFilterTask(unittest.TestCase):
    def setUp(self):
        pass

    def test_foo(self):
        self.assertEqual(1, 1)


    def get_valid_new_release_sequence_tsv(self):
        val = """PDB_ID  Sequence_Count  Sequence
                 2N1I    1       GSGFPTSEDFTPKEGSPYEAPVYIPEDIPIPADFELRESSIPGAGLGVWAKRKMEAGERLGPCVVVPRAAAKETDFGWEQILTDVEVSPQEGCITKISEDLGSEKFCVDANQAGAGSWLKYIRVACSCDDQNLTMCQISEQIYYKVIKDIEPGEELLVHVKEGVYPLGTVPPGLDE
                 2N27    1       ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK
                 2N3M    1       (DT)(DG)(DG)(DT)(DG)(DG)(DT)(DG)(DG)(DT)(DT)(DG)(DT)(DT)(DG)(DT)(DG)(DG)(DT)(DG)(DG)(DT)(DG)(DG)(DT)(DG)(DG)(DT)
                 2N4M    1       (DG)(DT)(DG)(DC)(4E9)(DT)(DG)(DT)(DT)(DT)(DG)(DT)
                 2N4M    2       (DA)(DC)(DA)(DA)(DA)(DC)(DA)(DC)(DG)(DC)(DA)(DC)
                 2N7A    1       MEGDRQYGDGYLLQVQELVTVQEGLSVHVPCSFSYPQDGWTDSDPVHGYWFRAGDRPYQDAPVATNNPDREVQAETQGRFQLLGDIWSNDCSLSIRDARKRDKGSYFFRLERGSMKWSYKSQLNYKTKQLSVFVTALTHGSLVPR
                 2N7B    1       MEGDRQYGDGYLLQVQELVTVQEGLSVHVPCSFSYPQDGWTDSDPVHGYWFRAGDRPYQDAPVATNNPDREVQAETQGRFQLLGDIWSNDCSLSIRDARKRDKGSYFFRLERGSMKWSYKSQLNYKTKQLSVFVTALTHGSLVPR
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

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()



