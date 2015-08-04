__author__ = 'robswift'

import os
import gzip
import socket
from d3r.PreRelease import filtering_sets as filtering_sets
from SequenceBase import Base
from collections import defaultdict
from Bio.Alphabet import generic_protein
from Bio import SeqIO


class Test(Base):
    """

    """
    if socket.gethostname() == 'Robs-MacBook-Pro.local':
        pdb_dir = '/Users/robswift/Documents/Work/D3R/devel/data/pdb'
    else:
        pdb_dir = '/data/pdb'

    pdb_dict = defaultdict(list)

    @staticmethod
    def set_pdb_dict(fasta):
        """
        Each PDB ID is mapped to a list of sequences, one sequence for each of its chains. This information is stored
        in a class-wide default dictionary pdb_dict, with the following structure
        pdb_dict = { 'PDB ID' : [ seq_chainA, seq_chainB, ... ] }
        :param fasta: path to the PDB sequences stored in FASTA format, i.e. "pdb_seqres.txt"
        :return:
        """
        fasta_handle = open(fasta, 'r')
        for record in SeqIO.parse(fasta_handle, "fasta"):
            pdb_id, chain_id = record.id.split('_')
            # only add protein sequences
            if 'mol:protein' in record.description:
                record.alphabet = generic_protein
                Test.pdb_dict[pdb_id].append(record)
        fasta_handle.close()

    def __init__(self):
        super(Base, self).__init__()
        self.sequences = {}
        self.coverage = None
        self.identity = None
        self.resolution = None
        self.ligands = { }                 # {'resname': ('InChI', 'label'}. label is 'do_not_call' or 'dock'
        self.exp_method = None

    def set_coverage(self, coverage):
        self.coverage = coverage

    def set_identity(self, identity):
        self.identity = identity

    def set_sequence(self, pdb_id, chain_id):
        """
        Looks for pdb sequence in the contents of the 'pdb_seqres.txt' file, stored in pdb_dict
        :param pdb_id: 4-letter pdb id
        :param chain_id: chain-id
        :return:
        """
        pdb_id = pdb_id.lower()  # PDBIDS are lower case in fasta
        self.sequences[pdb_id] = Test.pdb_dict[pdb_id]

    def set_resolution(self):
        """
        The resolution of the PDB of the target object is extracted from the PDB file stored in one of the
        subdirectories of Test.pdb_dir. If a PDB file cannot be found, a warning message is printed. If the resolution
        can not be found in the PDB file, a warning message is printed.
        """

        f = self.read_pdb_file()
        if f and self.exp_method == 'X-RAY DIFFRACTION':
            try:
                self.resolution = [x.split()[3] for x in f if 'REMARK   2 RESOLUTION' in x][0]
            except IndexError:
                print "Resolution for PDB %s could not be found" % self.pdb_id
        elif f and self.exp_method != 'X-RAY DIFFRACTION':
            self.resolution = 'Not applicable'
        else:
            print "PDB file %s could not be found. Resolution could not be set" % self.pdb_id

    def set_expt_method(self):
        """

        :return:
        """
        f = self.read_pdb_file()
        if f:
            try:
                self.exp_method = ' '.join([x.split()[1:] for x in f if 'EXPDTA' in x][0])
            except IndexError:
                print "The experimental method for PDB %s could not be found" % self.pdb_id
        else:
            print "PDB file %s could not be found. The experimental method could not be set" % self.pdb_id


    def set_ligands(self):
        """
        Add a ligand object to the ligand list (add a check to ensure it's a ligand class being passed)
        :param ligand:
        :return:
        """
        f = self.read_pdb_file()
        if f:
            for resname in [x.split()[1] for x in f if 'HETNAM' in x]:
                if not resname.isdigit():
                    if resname in filtering_sets.do_not_call:
                        label = 'do_not_call'
                    else:
                        label = 'dock'
                    self.ligands[resname] = ('NA', label)
        else:
            print "pdb file: %s could not be found. Test ligand dictionary could not be set" % self.pdb_id

    def read_pdb_file(self):
        """
        The PDB file corresponding to the PDB ID of the target object is read. Lines of the file are extracted into
        a list, f, and returned. If a PDB file cannot be found an empty list is returned.
        :return: f, a list of lines contained in the PDB file corresponding to the PDB ID of the target object.
        """
        pdb_file = os.path.join(Test.pdb_dir, self.pdb_id[1:3], 'pdb' + self.pdb_id + '.ent.gz')
        try:
            handle = gzip.open(pdb_file, 'rb')
            f = handle.readlines()
            handle.close()
        except IOError:
            f = []
        return f

    def get_coverage(self):
        return self.coverage

    def get_identity(self):
        return self.identity

    def get_sequences(self):
        return self.sequences

    def get_resolution(self):
        return self.resolution

    def get_ligand_names(self):
        return self.ligands.keys()

    def get_expt_method(self):
        return self.exp_method
