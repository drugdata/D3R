__author__ = 'robswift'

import os
import socket
from Bio.Alphabet import generic_protein
from Bio import SeqIO
from Bio.PDB import *
from Base import Base
from Ligand import Ligand
from d3r.filter import FilteringSets as filtering_sets
from d3r.blast.TestSequence import TestSequence

class Test(Base):
    """

    """
    if socket.gethostname() == 'Robs-MacBook-Pro.local':
        pdb_dir = '/Users/robswift/Documents/Work/D3R/devel/data/pdb'
    else:
        pdb_dir = '/data/pdb.extracted'

    pdb_dict = {}

    @staticmethod
    def set_pdb_dict(fasta):
        """
        Each PDB ID is mapped to a list of sequences, one sequence for each of its chains. This information is stored
        in a class-wide default dictionary pdb_dict, with the following structure
        pdb_dict = { 'pdbid_chainid' : Bio.SeqRecord }
        :param fasta: path to the PDB sequences stored in FASTA format, i.e. "pdb_seqres.txt"
        :return:
        """
        fasta_handle = open(fasta, 'r')
        for record in SeqIO.parse(fasta_handle, "fasta"):
            # only add protein sequences
            if 'mol:protein' in record.description:
                record.alphabet = generic_protein
                Test.pdb_dict[record.id] = record
        fasta_handle.close()

    def __init__(self):
        super(Test, self).__init__()
        self.retain = None                      # True if the test object is retained for docking, False if not.
        self.retain_reason = None               # A string that stores the reason the test object was retained for -->
                                                # --> docking.
        self.pdb = None                         # Bio.PDB structure object
        self.resolution = None                  # structural resolution (angstroms)
        self.exp_method = None
        self.highest_resolution = False         # Is this the highest resolution of a Target.test_list?
        self.most_similar = False               # Is one of the dockable ligands most similar to the Target ligand?
        self.least_similar = False              # Is one of the dockable ligands least similar to the Target ligand?

        self.retain_reasons = {
            1 : 'The BLAST hit is bound to the ligand with the largest maximum common substructure',
            2 : 'The BLAST hit is bound to the ligand with the smallest maximum common substructure',
            3 : 'The BLAST hit is the highest resolution structure',
            4 : 'The BLAST hit is the highest resolution apo structure',
        }

    def set_retain_reason(self, selection):
        """
        Sets the value of 'retain_reason', an attribute that explains why the test object is being retained for docking.
        In addition to setting a reason, this method also sets the value of 'retain' to True. To set a reason, an
        integer argument between 1 and 4, whose value maps to a particular reason, is provided. Selection to reason
        mapping:
            1 -> The BLAST hit is bound to the ligand with the largest maximum common substructure
            2 -> The BLAST hit is bound to the ligand with the smallest maximum common substructure
            3 -> The BLAST hit is the highest resolution structure
            4 -> The BLAST hit is the highest resolution apo structure
        Note that if the selection is outside of the range 1 to 4, no reason is set. Once the value for reason has
        been initiated, it is not fixed, and can be reset.
        :param selection: (int)
        """
        if int(selection) > 4 or int(selection) < 1:
            pass
        else:
            self.retain = True
            self.retain_reason = self.retain_reasons[int(selection)]

    def read_pdb(self):
        """
        Reads the PDB file corresponding to the Test objects wwPDB ID, self.pdb_id, and stores the results as a
        Bio.PDB structure object. If reading the file and assigning the Bio.PDB structure object is successful, True
        is returned. Otherwise, False is returned. Prior to calling the method, the pdb_id attribute must be set. E.g.
        typical usage follows
            test = d3r.Blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
        :return: Boolean
        """
        if not self.pdb_id:
            print "Test.pdb_id must be set prior to calling the read_pdb method"
            return False
        pdb_file = os.path.join(Test.pdb_dir, self.pdb_id[1:3], 'pdb' + self.pdb_id + '.ent')
        try:
            handle = open(pdb_file, 'r')
            handle.close()
        except IOError, e:
            print "Problems when assigning the ligand for wwPDB ID: %s\n%s" % (self.pdb_id, e)
            return False
        parser = PDBParser()
        try:
            self.pdb = parser.get_structure(format(self.pdb_id), pdb_file)
            if not self.pdb:
                return False
            else:
                return True
        except:
            return False

    def set_coverage(self, coverage):
        self.coverage = float(coverage)

    def set_identity(self, identity):
        self.identity = float(identity)

    def set_sequence(self, pdb_id, chain_id, alignment, query_length):
        """
        Looks for pdb sequence in the contents of the 'pdb_seqres.txt' file, stored in pdb_dict
        :param pdb_id: 4-letter pdb id
        """
        id = pdb_id.lower() + '_' + chain_id.upper()
        seq = TestSequence(Test.pdb_dict[id])
        seq.pdb_id = pdb_id
        seq.chain_id = chain_id
        seq.query_length = query_length
        seq.alignment = alignment
        seq.set_coverage_and_identity()
        self.sequences.append(seq)
        self.chain_count = len(self.sequences)

    def set_ligands(self):
        """
        Creates a d3r.blast.Ligand object for each of the bound ligands found in the PDB, and appends
        this object to the appropriate list, either dock or do_not_call. If the ligand object is added
        to the dock list, the dock_count and the ligand_count is incremented. If the ligand object is added to the
        do_not_call list, then the ligand_count is incremented. No effort is made to check whether or not the
        input ligand is unique, i.e. that it hasn't been added before. Prior to calling the method, the pdb_id
        attribute must be set, and the read_pdb method must be called. For example:
            test = d3r.Blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
            test.set_ligands()
        """
        res_list = Selection.unfold_entities(self.pdb, 'R')
        hetero_list = [res for res in res_list if 'H_' in res.id[0] and res.id[0] != 'H_MSE']
        for hetero in hetero_list:
            ligand = Ligand()
            ligand.resname = hetero.resname
            if hetero.resname in filtering_sets.do_not_call:
                ligand.label = 'do_not_call'
                self.do_not_call.append(ligand)
            else:
                ligand.label = 'dock'
                ligand.set_rd_mol_from_bio_pdb_res(hetero)
                self.dock.append(ligand)
                self.dock_membership[ligand.resname] = len(self.dock) - 1

    def set_resolution(self):
        """
        Set the resolution of the PDB of the test object, self.resolution. Before this method is called, the pdb_id
        attribute must be set, and the read_pdb method must be called. For example:
            test = d3r.Blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
            test.set_resolution()
        If the resolution cannot be set, False is returned, otherwise True is returned.
        :return: Boolean
        """
        if not self.pdb_id or not self.pdb:
            print "Before the set_resolution method is called, the pdb_id attribute must be set, and the read_pdb " \
                  "method must be called."
            return False
        # TO DO: add error checking for NMR structures, and other experimental methods from which resolution is
        # unavailable
        try:
            self.resolution = self.pdb.header['resolution']
            if self.resolution:
                return True
            else:
                return False
        except:
            return False

    def set_expt_method(self):
        """
        Set the experimental structure determination method of the PDB of the test object, self.exp_method. Before this
        method is called, the pdb_id attribute must be set, and the read_pdb method must be called. For example:
            test = d3r.Blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
            test.set_expt_method()
        If the resolution cannot be set, False is returned, otherwise True is returned.
        :return: Boolean
        """
        if not self.pdb_id or not self.pdb:
            print "Before the set_expt_method method is called, the pdb_id attribute must be set, and the read_pdb " \
                  "method must be called."
            return False
        try:
            self.exp_method = self.pdb.header['structure_method']
            if self.resolution:
                return True
            else:
                return False
        except:
            return False

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
