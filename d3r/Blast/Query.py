__author__ = 'robswift'

import os
from StringIO import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from rdkit import Chem
from rdkit.Chem import rdFMCS
from Base import Base
from Hit import Hit
from Ligand import Ligand

class Query(Base):
    """
    This class contains information about the macromolecular structures uniquely identified by a wwPDB in the
    new_release_structure_sequence.tsv and new_release_structure_nonpolymer.tsv files. The Target class contains
    a list of d3r.Test objects that are instantiated following a BLASTP of the Target sequence.
    Information is included to facilitate filtering:
      (i) self.contains_largest_mcss: Set to True if the maximum common substructure between the dockable ligand of the
      Target object, self.dock, and one of the contained d3r.blast.Test objects dockable ligands is largest among a
      set of Target objects.
      (ii) self.contains_smallest_mcss: Set to True if the maximum common substructure between the dockable ligand of
      the Target object, self.dock, and one of the contained d3r.blast.Test objects dockable ligands is smallest among
      a set of Target objects.
      (iii) self.inchi_error: Set to True if the Target object contains a dockable ligand and an error is produced when
      attempting to convert the InChI to a rdkit Mol object. Set to False if the Target object contains a dockable
      ligand and the InChI was successfully converted. Set to None if the Target object is Apo, or contains no dockable
      ligand.
    """
    def __init__(self):
        super(Query, self).__init__()             # inherit superclass constructor
        self.test_list = []                        # [test object, test object, ... ]
        self.test_membership = {}                  # {'PDB_ID' : index (int)} where PDB_ID is a lowercase string and
                                                   # index is the index of the corresponding PDB in test_list
        self.inchi_error = None                    # None, True, or False.
        self.contains_largest_mcss = None          #

    def set_ligand(self, resname, inchi, label):
        """
        Creates a d3r.blast.Ligand object entry for the ligand. If the ligand is dockable (not in in the
        d3r.filter.FilteringSets do_not_call set), the Ligand object is appended to the dock list, and the dock_count
        and ligand_count attributes are each incremented. If ligand is in the do_not_call filtering set, then the
        ligand is added to the do_not_call list, and the ligand_count attribute is incremented. Note that, no effort is
        made to check if the input ligand is unique or has not been added before. In the case of a dockable ligand, an
        attempt is made to convert the InChI to a rd mol object. If the conversion is successful, then the Target
        inchi_error attribute is set to False. If the conversion is unsuccessful, then the the inchi_error attribute is
        set to True. Non-dockable ligands are not converted to rd mol objects.
        :param resname: (string) the PDB resname of the input ligand
        :param inchi: (string) the inchi string, which uniquely identifies the ligand
        :param label: (string) set to either 'dock' or 'do_not_call'
        """
        ligand = Ligand(resname, inchi)
        if label == 'dock':
            if ligand.set_rd_mol_from_inchi():
                self.inchi_error = False
            else:
                self.inchi_error = True
            self.dock.append(ligand)
            self.dock_count += 1
            self.ligand_count += 1
        elif label == 'do_not_call':
            self.do_not_call.append(ligand)
            self.ligand_count += 1
        else:
            pass

    def set_sequence(self, chain_id, seq):
        """
        Converts an input FASTA sequence to a Bio.SeqRecord with the id field set to the input integer chain id. Each
        unique sequence receives a unique integer chain id, so chain ids will be numbered 1 through N, where N is the
        count of the last unique chain of the protein. Additionally, each time a sequence is added, the
        chain_count attribute of the Target object will be incremented. Note that no effort is made to check whether or
        not the input sequence is unique, i.e. that it hasn't been added before.
        :param chain_id: (an integer cast
        :param seq:  a FASTA sequence
        """
        self.sequences.append(SeqRecord(Seq(seq, generic_protein), id = format(chain_id)))
        self.chain_count += 1

    def run_blast(self, pdb_db, out_dir):
        """
        Runs BLASTP on the sequence(s) contained in the Target instance. If the chain of a sequence has not been
        assigned to the Target instance, then the sequence returns False. Otherwise, the sequence returns a
        Bio.Blast.Record object. If the Target instance is a multimer, the alignment objects in the list,
        Bio.Blast.Record.alignments are the intersection of the alignment objects, or BLASTP hits that result from
        running BLASTP on each of the sequences of each of the unique chains that constitutes the multimer. If the
        Target instance is a monomer, the alignment objects in the list Bio.Blast.Record.alignments is simply the list
        of BLASTP hits that result from running BLASTP on the single unique sequence that makes up the monomer.
        :param pdb_db: (string) The absolute path to a BLASTP database
        :param out_dir: (string) The absolute path to the output directory. A fasta file is writen here prior to running
        BLASTP and is deleted after the BLASTP search is completed.
        :return: Bio.Blast.Record object
        """
        if self.chain_count == 0:
            print "Fasta sequences have not been added to the Target.sequences list, and BLAST cannot be run"
            return False
        elif self.chain_count == 1:
            record = self.blast_monomer(pdb_db, out_dir)
            return record
        elif self.chain_count > 1:
            record = self.blast_multimer(pdb_db, out_dir)
            return record

    def blast_monomer(self, pdb_db, out_dir):
        """
        Runs a BLASTP search of sequences in the wwPDB for a single monomer chain
        :param pdb_db: (string) The absolute path to a BLASTP database
        :param out_dir: (string) The absolute path to the output directory. A fasta file is writen here prior to running
        :return: Bio.Blast.Record
        """
        for sequence in self.sequences:
            fasta = self.write_fasta(sequence, out_dir)
            if fasta:
                record = self.blastp(fasta, pdb_db)
                os.remove(fasta)
                return record

    def blast_multimer(self, pdb_db, out_dir):
        """
        Runs a BLASTP search of sequences in the wwPDB for each unique chain in a multimer.
        :param pdb_db: (string) The absolute path to a BLASTP database
        :param out_dir: (string) The absolute path to the output directory. A fasta file is writen here prior to running
        :return: Bio.Blast.Record
        """
        records = []
        for sequence in self.sequences:
            fasta = self.write_fasta(sequence, out_dir)
            if fasta:
                records.append(self.blastp(fasta, pdb_db))
                os.remove(fasta)
        if records:
            record = self.get_intersection(records)
            return record

    def get_intersection(self, records):
        """
        Returns a Bio.Blast record object, whose alignment objects are the intersection of the alignments stored in
        each of the Bio.Blast record objects contained in the input records list.
        :param records: a list of Bio.Blast.Record objects
        :return: record, a Bio.Blast.Record object
        """
        ids=[]
        for record in records:
            ids.append(set([self.get_pdb_id_from_alignment(a) for a in record.alignments]))
        id_intersection = set.intersection(*ids)
        self.clean_alignments(records[0], id_intersection)
        return records[0]

    def clean_alignments(self, record, id_intersection):
        """
        Removes all alignments from the input Bio.Blast.Record object whose sequences correspond to wwPDB IDs not in the
        input intersection set. It's assumed that the wwPDB IDs of the input Bio.Blast.Record.alignments list can form
        a superset of the input intersection set.
        :param record: a Bio.Blast.Record object
        :param id_intersection: (set) a set of wwPDB ids that represent the intersection of two or more BLASTP results.
        :return:
        """
        delete = []
        for index in range(len(record.alignments)):
            if self.get_pdb_id_from_alignment(record.alignments[index]) not in id_intersection:
                delete.append(index)
        delete.sort(reverse=True)
        for index in delete:
            del record.alignments[index]

    def blastp(self, fasta, pdb_db):
        """
        Runs BLASTP locally on a input fasta file and specified BLASTP database
        :param fasta: The absolute path to a FASTA file
        :param pdb_db: A BLASTP database
        :return: Bio.Blast.Record object
        """
        cline = NcbiblastpCommandline(cmd='blastp', query=fasta, db=pdb_db, evalue=0.001, outfmt=5)
        std_out, std_err = cline()
        blast_records = NCBIXML.parse(StringIO(std_out))
        record = next(blast_records)
        return record

    def write_fasta(self, sequence, out_dir):
        """
        Writes the input BioPython sequence records to a fasta file in the specified directory. If the fasta file was
        written successfully, the absolute path to the fasta file is returned. Otherwise False is returned.
        :sequence: (Bio.SeqRecord)
        :param outdir: (string) The absolute path of the directory where the fasta file will be written
        return: the absolute path to the fasta file
        """
        fasta = os.path.join(os.path.abspath(out_dir), '_'.join((self.pdb_id, sequence.id)) + '.fasta')
        try:
            SeqIO.write(sequence, fasta, "fasta")
            return fasta
        except:
            return False

    def get_pdb_id_from_alignment(self, alignment):
        """
        Returns a string of the four letter PDB ID found in alignment
        :param alignment:
        :return: 'pdb id'
        """
        pdb_chain_id = alignment.hit_def.encode('ascii').split()[0]
        pdb_id = pdb_chain_id.split('_')[0].lower()
        return pdb_id

    def get_chain_id_from_alignment(self, alignment):
        """
        Returns a string of the chain id, e.g. 'A', 'B', etc.
        :param alignment:
        :return:
        """
        pdb_chain_id = alignment.hit_def.encode('ascii').split()[0]
        chain = pdb_chain_id.split('_')[1].upper()
        return chain

    def sort_by_resolution(self):
        """
        If a test list exists, and each test object has a resolution value, the list is sorted from highest resolution
        (lowest value) to lowest resolution (highest value)
        """
        self.test_list.sort(key=lambda x: x.get_resolution())

    def set_test_list(self, record):
        """
        Populate self.test_list with a list of test objects. Each test object corresponds to a BLASTP hit stored in one
        of the alignment objects found in the input record.alignments list and contains the following information about
        the structural hit: wwPDB id, % coverage (alignment length / query length), % identity (number of identical
        amino acids in the alignment / hit length), experimental structural determination method, resolution,
        :param record (Bio.Blast.Record object)
        """
        for alignment in record.alignments:
            pdb_id = self.get_pdb_id_from_alignment(alignment)
            chain_id = self.get_chain_id_from_alignment(alignment)
            if pdb_id in self.test_membership.keys():
                self.test_list[self.test_membership[pdb_id]].set_sequence(pdb_id, chain_id, alignment,
                                                                          record.query_length)
            else:
                test = Hit()
                # set the pdb_id
                test.pdb_id = pdb_id
                # read in a pdb file into the test object
                test.read_pdb()
                # set the resolution
                test.set_resolution()
                # set the experimental method
                test.set_expt_method()
                # set the test sequence
                test.set_sequence(pdb_id, chain_id, alignment, record.query_length)
                # add ligands
                test.set_ligands()
                # update the test_list
                self.test_list.append(test)
                # update test_membership
                self.test_membership[pdb_id] = len(self.test_list) - 1


    def mcss(self):
        """
        Performs a maximum common substructure (MCSS) search between each of the Target objects dockable ligands and
        each of the Test objects dockable ligands. If the calculation is successful, the MCSS is stored as an rd mol
        object in the Ligand embedded in the Test object.
        each of the dockable ligands of each of the corresponding Test structures, i.e. those in self.test_list.
        """
        if self.dock_count > 0:
            for target_lig in self.dock:
                for test in self.test_list:
                    if test.dock_count > 0:
                        test_index = self.test_membership[test.pdb_id]
                        for test_lig in test.dock:
                            dock_index = test.dock_membership[test_lig.resname]
                            res = rdFMCS.FindMCS([target_lig, test_lig])
                            if res.numAtoms == 0:
                                self.test_list[test_index].dock[dock_index].mcss_error = True
                            else:
                                mcss = Chem.MolFromSmarts(res.smartsString)
                                self.test_list[test_index].dock[dock_index].mcss = mcss

