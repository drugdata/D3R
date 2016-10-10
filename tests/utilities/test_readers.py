__author__ = 'j-wags'

import unittest
import tempfile
import os.path
import shutil 
from d3r.utilities.readers import ReadText


"""
test_readers
--------------------------------

Tests for the 'readers' module
"""

class TestReaders(unittest.TestCase):
    def setUp(self):
        self.test_txt_filename = self.make_test_txt()

    def make_test_txt(self, output_filename = '5aqc.txt'):
        data = '''query, 5aqc
ph, 7
ligand, 5JB
inchi, InChI=1S/C48H76N7O18P3S/c1-27(32-12-13-33-31-11-10-29-22-30(56)14-17-47(29,5)34(31)15-18-48(32,33)6)8-7-9-28(2)45(61)77-21-20-50-36(57)16-19-51-43(60)40(59)46(3,4)24-70-76(67,68)73-75(65,66)69-23-35-39(72-74(62,63)64)38(58)44(71-35)55-26-54-37-41(49)52-25-53-42(37)55/h22,25-28,31-35,38-40,44,58-59H,7-21,23-24H2,1-6H3,(H,50,57)(H,51,60)(H,65,66)(H,67,68)(H2,49,52,53)(H2,62,63,64)/t27-,28+,31+,32-,33+,34+,35-,38-,39-,40+,44-,47+,48-/m1/s1
size, 88
rotatable_bond, 25
LMCSS, 5cxi, 5TW, chain: A, (size: 82, mcss_size: 81, resolution:  2.0) 
SMCSS, 5cw8, 55X, chain: A, (size: 84, mcss_size: 58, resolution:  2.6) 
hiResHolo, 5cxi, 5TW, chain: A, (resolution:  2.0)
hiResApo, 3mnl
hiTanimoto, 5cxi, 5TW, chain: A, (tanimoto_similarity: 0.98, resolution:  2.0)
hiTanimoto, 5cw8, 55X, chain: A, (tanimoto_similarity: 0.97, resolution:  2.6)
'''
        with open(output_filename, 'wb') as of:
            of.write(data)
        return output_filename


    def test_read_text(self):
        test_filename = self.make_test_txt()
        readtext_obj = ReadText()
        readtext_result = readtext_obj.parse_txt(test_filename)
        known_answer = self.get_parsed_test_data()
        for key in known_answer.keys():
            self.assertEqual(readtext_result[key], known_answer[key])



    def get_parsed_test_data(self):
        parsed_dict = {'ligand': ['5JB'], 
                       'rotatable_bond': ['25'], 
                       'size': ['88'],
                       'inchi': ['InChI=1S/C48H76N7O18P3S/c1-27(32-12-13-33-31-11-10-29-22-30(56)14-17-47(29,5)34(31)15-18-48(32,33)6)8-7-9-28(2)45(61)77-21-20-50-36(57)16-19-51-43(60)40(59)46(3,4)24-70-76(67,68)73-75(65,66)69-23-35-39(72-74(62,63)64)38(58)44(71-35)55-26-54-37-41(49)52-25-53-42(37)55/h22,25-28,31-35,38-40,44,58-59H,7-21,23-24H2,1-6H3,(H,50,57)(H,51,60)(H,65,66)(H,67,68)(H2,49,52,53)(H2,62,63,64)/t27-,28+,31+,32-,33+,34+,35-,38-,39-,40+,44-,47+,48-/m1/s1'], 
                       'query': ['5aqc'], 
                       'ph': ['7'], 
                       'LMCSS': [{'lig_name': ['5TW'], 
                                  'chain': ['A'], 
                                  'cand_id': ['5cxi'], 
                                  'mcss_size': ['81'], 
                                  'resolution': ['2.0'], 
                                  'size': ['82']}], 
                       'hiResApo': [{'cand_id': ['3mnl']}], 
                       'hiResHolo': [{'cand_id': ['5cxi'], 
                                      'lig_name': ['5TW'], 
                                      'resolution': ['2.0'], 
                                      'chain': ['A']}], 
                       'SMCSS': [{'lig_name': ['55X'], 
                                  'chain': ['A'], 
                                  'cand_id': ['5cw8'], 
                                  'mcss_size': ['58'], 
                                  'resolution': ['2.6'], 
                                  'size': ['84']}], 
                       'hiTanimoto': [{'cand_id': ['5cxi'], 
                                       'lig_name': ['5TW'], 
                                       'resolution': ['2.0'], 
                                       'tanimoto_similarity': ['0.98'], 
                                       'chain': ['A']}, 
                                      {'cand_id': ['5cw8'], 
                                       'lig_name': ['55X'], 
                                       'resolution': ['2.6'], 
                                       'tanimoto_similarity': ['0.97'], 
                                       'chain': ['A']}
                                      ], 
                       }

        return parsed_dict

    def tearDown(self):
        os.remove(self.test_txt_filename)


if __name__ == '__main__':
    unittest.main()
