__author__ = 'robswift'

from collections import defaultdict

class Base(object):
    def __init__(self):
        self.pdb_id = None                 # 4 letter PDB code
        self.sequences = defaultdict(list) # {'id' : [SeqRecord, SeqRecord] }. Where chain_id is 'pdb_id_chain_id'

    # setters
    def set_pdb_id(self, pdb_id):
        self.pdb_id = pdb_id.lower()

    # getters
    def get_pdb_id(self):
        return self.pdb_id
