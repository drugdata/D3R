__author__ = 'robswift'

import os
import sys
import re
import logging

from Bio.Alphabet import IUPAC
from Bio import SeqIO
from d3r.blast.base import Base
from d3r.blast.ligand import Ligand
from d3r.filter import filtering_sets as filtering_sets
from d3r.blast.hit_sequence import HitSequence

logger = logging.getLogger(__name__)

try:
    from Bio.PDB import *
except ImportError:
    logger.exception('Unable to import Bio.PDB Hit class may not work')


class RegDict(dict):
    def get_matching(self, event):
        return ((key, self[key]) for key in self if re.match(event, key))


class Hit(Base):
    """

    """

    pdb_dir = None
    pdb_dict = RegDict()

    @staticmethod
    def set_pdb_dir(pdb_path):
        if os.path.isdir(pdb_path):
            Hit.pdb_dir = pdb_path
        else:
            raise IOError("'{pdb_path}' is not a directory".format(pdb_path=pdb_path))

    @staticmethod
    def set_pdb_dict(fasta):
        """
        Each PDB ID is mapped to a list of sequences, one sequence for each of its chains. This information is stored
        in a class-wide default dictionary pdb_dict, with the following structure
        pdb_dict = { 'pdbid_chainid' : Bio.SeqRecord }
        :param fasta: path to the PDB sequences stored in FASTA format, i.e. "pdb_seqres.txt"
        """
        try:
            fasta_handle = open(fasta, 'r')
        except IOError:
            logger.exception('Caught exeption opening fasta file')
            sys.exit(1)
        logger.debug('Parsing ' + fasta)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            # only add protein sequences
            if 'mol:protein' in record.description:
                record.seq.alphabet = IUPAC.protein
                Hit.pdb_dict[record.id] = record
            elif 'mol:na' in record.description:
                Hit.pdb_dict[record.id] = record
                # default SingleLetterAlphabet() assigned

        fasta_handle.close()

    def __init__(self):
        super(Hit, self).__init__()
        self.retain = None              # True if the test object is retained for docking, False if not.
        self.reasons_to_retain = []     # A list that stores the reasons the test object was retained for docking.
        self.pdb = None                 # Bio.PDB structure object
        self.chain_count = 0            # an integer count of the number of PDB chains, which may correspond to -->
                                        # --> identical sequences.
        self.resolution = None          # structural resolution (angstroms)
        self.exp_method = None          # Method used to determine the structure, e.g. x-ray diffraction'
        self.largest_mcss = None        # The the largest MCSS contained in the self.dock.mcsss list
        self.largest_mcss_chain = []    # The chain ID of the largest MCSS
        self.largest_index = []         # The index of the largest MCSS in the self.dock list
        self.smallest_mcss = None       # The smallest MCSS contained in the self.dock.mcsss list
        self.smallest_mcss_chain = []   # The chain ID of the smallest MCSS
        self.smallest_index = []        # The index of the smallest MCSS in the self.dock list
        self.highest_tanimoto = None    # The tanimoto similarity contained in the self.dock.mcsss list
        self.highest_tanimoto_chain = []# The chain ID of the highest tanimoto
        self.highest_tanimoto_index = []# The index of the highest tanimoto in the self.dock list
        self.sequence_membership = {}   # {'pdb_id'_'chain_id' : index} where index is the index of the -->
                                        # --> corresponding HitSequence object in the sequences list.

        #add category 5 for largest tanimoto score 09/01 sliu. 
        self.retain_reasons = {
            1: 'The BLAST hit is bound to the ligand with the largest maximum common substructure',
            2: 'The BLAST hit is bound to the ligand with the smallest maximum common substructure',
            3: 'The BLAST hit is the highest resolution holo structure',
            4: 'The BLAST hit is the highest resolution apo structure',
            5: 'The BLAST hit is bound to the ligand with the highest tanimoto score in structural similarity',
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
            5 -> The BLAST hit is bound to the ligand with the highest tanimoto score in structural similarity
        Note that if the selection is outside of the range 1 to 5, no reason is set. Once the value for reason has
        been initiated, it is not fixed, and can be reset.
        :param selection: (int)
        """
        #add category 5 for largest tanimoto score 09/01 sliu. 
        if int(selection) > 5 or int(selection) < 1:
            pass
        else:
            self.retain = True
            reason = self.retain_reasons[int(selection)]
            if reason not in self.reasons_to_retain:
                self.reasons_to_retain.append(reason)

    def set_maxmin_mcss(self):
        """
        Finds the largest and smallest MCSS objects contained in the list of dockable ligands. If this list of ligands
        doesn't exist, or if none of the ligands contain an MCSS, the method does not execute.
        """
        if not self.dock:
            pass
        elif not self.dock[0].mcsss:
            pass
        else:
            self.largest_mcss = self.dock[0].mcsss[0]
            self.smallest_mcss = self.dock[0].mcsss[0]
            #add highest tanimoto score 09/01 sliu
            self.highest_tanimoto = self.dock[0].mcsss[0]
            self.largest_mcss_chain = []
            self.smallest_mcss_chain = []
            self.highest_tanimoto_chain = []
            #store the index of the largest mcss in the self.dock object 
            self.largest_index = []
            #store the index of the smallest mcss in the self.dock object 
            self.smallest_index = [] 
            self.highest_tanimoto_index = [] 
            for (ligand_index, ligand) in enumerate (self.dock):
                    #modified by sliu 08/04 to add chain info
                    for mcss in ligand.mcsss:
                        if mcss.size > self.largest_mcss.size:
                            self.largest_mcss = mcss
                            #if this mcss is bigger then clean up the chain list and pick this one  
                            self.largest_mcss_chain = []
                            self.largest_mcss_chain.append(ligand.chain)
                            self.largest_index = []
                            self.largest_index.append(ligand_index)
                        if mcss.size == self.largest_mcss.size:
                            #if this mcss is the same as the original one, just append
                            self.largest_mcss_chain.append(ligand.chain)
                            self.largest_index.append(ligand_index)
                        if mcss.size < self.largest_mcss.size:
                            self.smallest_mcss = mcss
                            self.smallest_mcss_chain = []
                            self.smallest_mcss_chain.append(ligand.chain)      
                            self.smallest_index = []                                
                            self.smallest_index.append(ligand_index)
                        if mcss.size == self.smallest_mcss.size:
                            self.smallest_mcss_chain.append(ligand.chain)
                            self.smallest_index.append(ligand_index)
                        if mcss.tanimoto > self.highest_tanimoto.tanimoto:
                            self.highest_tanimoto = mcss
                            self.highest_tanimoto_chain = []
                            self.highest_tanimoto_chain.append(ligand.chain)
                            self.highest_tanimoto_index = []
                            self.highest_tanimoto_index.append(ligand_index)
                        if mcss.tanimoto == self.highest_tanimoto.tanimoto:
                            self.highest_tanimoto_chain.append(ligand.chain)
                            self.highest_tanimoto_index.append(ligand_index)
                         
            #sort the chain list so that later we could anyway pick the first one.
            (self.largest_mcss_chain, self.largest_index) = (list(x) for x in zip(*sorted(zip(self.largest_mcss_chain, self.largest_index))))
            #print "CCCCCCCCCCCCCCC", (self.largest_mcss_chain, self.largest_index)  
            (self.smallest_mcss_chain, self.smallest_index) = (list(x) for x in zip(*sorted(zip(self.smallest_mcss_chain, self.smallest_index))))
            #here we need to clean up the self.dock to only reserve one ligand
            (self.highest_tanimoto_chain, self.highest_tanimoto_index) = (list(x) for x in zip(*sorted(zip(self.highest_tanimoto_chain, self.highest_tanimoto_index))))

    def read_pdb(self):
        """
        Reads the PDB file corresponding to the Test objects wwPDB ID, self.pdb_id, and stores the results as a
        Bio.PDB structure object. If reading the file and assigning the Bio.PDB structure object is successful, True
        is returned. Otherwise, False is returned. Prior to calling the method, the pdb_id attribute must be set. E.g.
        typical usage follows
            test = d3r.blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
        :return: Boolean
        """
        logger.debug('In read_pdb()')
        if not self.pdb_id:
            logger.error("Test.pdb_id must be set prior to calling the read_pdb method")
            return False
        pdb_file = os.path.join(Hit.pdb_dir, self.pdb_id[1:3], 'pdb' + self.pdb_id + '.ent')
        try:
            handle = open(pdb_file, 'r')
            handle.close()
        except IOError as e:
            logger.exception('Problems when assigning the ligrand for wwPDB ID: ' + self.pdb_id)
            return False
        try:
            parser = PDBParser()
            self.pdb = parser.get_structure(format(self.pdb_id), pdb_file)
            if not self.pdb:
                logger.error('Unable to get structure for pdb')
                return False
            else:
                return True
        except:
            logger.exception('Caught exception')
            return False

    def set_coverage(self, coverage):
        self.coverage = float(coverage)

    def set_identity(self, identity):
        self.identity = float(identity)

    def set_sequence(self, pdb_id, chain_id, record, alignment):
        """
        Looks for pdb sequence in the contents of the 'pdb_seqres.txt' file, stored in pdb_dict
        :param pdb_id: 4-letter pdb id
        """

        id = '{pdb_id}_{chain_id}'.format(pdb_id = pdb_id.lower(), chain_id = chain_id.upper())
        logger.debug('Looking for ' + id + ' in pdb_seqres.txt')
        if id in self.sequence_membership.keys():
            self.sequences[self.sequence_membership[id]].set_query_alignment(record, alignment)
        else:
            seq_record = Hit.pdb_dict[id]
            if self.unique_sequence(seq_record.seq):
                self.sequence_count += 1
                hs = HitSequence(seq_record)
                hs.hit_sequence_id = self.sequence_count
            else:
                seq_id = [hs.hit_sequence_id for hs in self.sequences if hs.seq_record.seq == seq_record.seq][0]
                hs = HitSequence(seq_record)
                hs.hit_sequence_id = seq_id
            hs.blast_hit = True
            hs.hit_pdb_id = pdb_id
            hs.hit_chain_id = chain_id
            hs.set_query_alignment(record, alignment)
            self.sequences.append(hs)
            self.sequence_membership[id] = len(self.sequences) - 1
            self.chain_count += 1

    def fill_sequence(self):
        """
        Fills in chains present in the pdb that weren't picked up by blast. For example, nucleic acid
        components co-crystallized with the protein would be picked up by the fill_sequence method. This should only
        be called after the set_hits method is called
        """
        logger.debug('In fill_sequence()')
        for id, seq_record in Hit.pdb_dict.get_matching(self.pdb_id):
            if id in self.sequence_membership.keys():
                # the chain exists in the sequence list, so no need to do anything
                continue
            else:
                # the chain hasn't been assigned. There are two possibilities here.
                # i) the chain corresponds to a sequence that has already been assigned.
                # ii) the chain corresponds to a new sequence
                # case i) is dealt with first
                if self.unique_sequence(seq_record.seq):
                    self.sequence_count += 1
                    hs = HitSequence(seq_record)
                    hs.hit_sequence_id = self.sequence_count
                else:
                    seq_id = [hs.hit_sequence_id for hs in self.sequences if hs.seq_record.seq == seq_record.seq][0]
                    hs = HitSequence(seq_record)
                    hs.hit_sequence_id = seq_id
                hs.blast_hit = False
                pdb_id, chain_id = id.split('_')
                hs.hit_pdb_id = pdb_id
                hs.hit_chain_id = chain_id
                self.sequences.append(hs)
                self.sequence_membership[id] = len(self.sequences) - 1
                self.chain_count += 1

    def unique_sequence(self, seq):
        """
        Returns true if the input sequence cannot be found in the list of Hit object's sequences, otherwise returns
        False.
        :param seq: a Bio.Seq object
        :return: Boolean
        """
        if len([hs for hs in self.sequences if hs.seq_record.seq == seq]) == 0:
            return True
        else:
            return False

    def set_ligands(self, chain_id):
        """
        Creates a d3r.blast.Ligand object for each of the bound ligands found in the PDB, and appends
        this object to the appropriate list, either dock or do_not_call. If the ligand object is added
        to the dock list, the dock_count and the ligand_count is incremented. If the ligand object is added to the
        do_not_call list, then the ligand_count is incremented. No effort is made to check whether or not the
        input ligand is unique, i.e. that it hasn't been added before. Prior to calling the method, the pdb_id
        attribute must be set, and the read_pdb method must be called. For example:
            test = d3r.blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
            test.set_ligands()
        """
        logger.debug('In set_ligands()')
        model_list = Selection.unfold_entities(self.pdb, 'M')
        chain_list = Selection.unfold_entities(model_list[0], 'C')
        chain_format = "<Chain id=%s>"%chain_id
        chain_index = 0 
        for (index, chain_object) in enumerate(chain_list):
            if str(chain_object) == chain_format:
                chain_index = index
        #print "Chain index", chain_index
        res_list = Selection.unfold_entities(chain_list[chain_index], 'R')
        hetero_list = [res for res in res_list if 'H_' in res.id[0] and res.id[0] != 'H_MSE']
        #modified by sliu 08/04, add chain info for ligand
        for hetero in hetero_list:
            assigned = []
            for l in self.dock:
                ligand_info = (l.resname, l.chain)
                assigned.append(ligand_info)
            for l in self.do_not_call:
                ligand_info = (l.resname, l.chain)
                assigned.append(ligand_info)

            resname = hetero.resname.strip()
            if (resname, chain_id) in assigned:
                #the resname, chain id info is already in stored
                continue
            ligand = Ligand()
            ligand.resname = resname
            ligand.chain = chain_id
            if resname in filtering_sets.do_not_call:
                ligand.label = 'do_not_call'
                self.do_not_call.append(ligand)
                self.ligand_count += 1
            else:
                ligand.label = 'dock'
                ligand.set_rd_mol_from_resname(resname)
                #modified by sliu 08/08
                ligand.set_size()
                ligand.set_heavy_size()
                self.dock.append(ligand)
                self.dock_count += 1
                self.ligand_count += 1

    def set_resolution(self):
        """
        Set the resolution of the PDB of the test object, self.resolution. Before this method is called, the pdb_id
        attribute must be set, and the read_pdb method must be called. For example:
            test = d3r.blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
            test.set_resolution()
        If the resolution cannot be set, False is returned, otherwise True is returned.
        :return: Boolean
        """
        logger.debug('In set_resolution()')
        if not self.pdb_id or not self.pdb:
            logger.error("Before the set_resolution method is called, the pdb_id attribute must be set, "
                         "and the read_pdb "
                         "method must be called.")
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
            logger.exception('Caught exception')
            return False

    def set_expt_method(self):
        """
        Set the experimental structure determination method of the PDB of the test object, self.exp_method. Before this
        method is called, the pdb_id attribute must be set, and the read_pdb method must be called. For example:
            test = d3r.blast.Test()
            test.pdb_id = '1xdn'
            test.read_pdb()
            test.set_expt_method()
        If the resolution cannot be set, False is returned, otherwise True is returned.
        :return: Boolean
        """
        logger.debug('In set_expt_method()')
        if not self.pdb_id or not self.pdb:
            logger.error("Before the set_expt_method method is called, the pdb_id attribute must be set, and the read_pdb "
                         "method must be called.")
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

