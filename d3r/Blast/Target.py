__author__ = 'robswift'

from collections import defaultdict
import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from StringIO import StringIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import ParseAlignment
from SequenceBase import Base
from Test import Test


class Target(Base):
    """

    """


    def __init__(self):
        super(Base, self).__init__()                # inherit superclass constructor
        self.sequences = defaultdict(list)          # {'pdbid' : [SeqRecord, ..., SeqRecord]}
        self.test_list = [ ]                        # {'pdbid' : [test objects]}
        self.ligands = { }                          # {'resname': ('InChI', 'label'}. label is 'do_not_call' or 'dock'

    def set_ligand(self, resname, inchi, label):
        """

        :param ligand:
        :return:
        """
        self.ligands[resname] = (inchi, label)

    def set_sequence(self, pdb_id, chain_id, seq):
        id = pdb_id + '_' + chain_id
        self.sequences[pdb_id].append(SeqRecord(Seq(seq, generic_protein), id = id))

    def get_ligand_names(self):
        """
        Return ligand dictionary
        :return:
        """
        return self.ligands.keys()

    def get_sequences(self):
        return self.sequences

    def run_blast(self, pdb_db, out_dir):
        """
        Run BLASTP
        :return:
        """
        # this will only work for the simplest monomer case....add functionality for oligomer
        for key in self.sequences.keys():
            # write a fasta file
            fasta_file = os.path.join(os.path.abspath(out_dir), key + '.fasta')
            SeqIO.write(self.sequences[key], fasta_file, "fasta")
            #FastaIO.write_fasta(fasta_file, self.sequences[key])

            # run blastp
            cline = NcbiblastpCommandline(cmd='blastp', query=fasta_file, db=pdb_db, evalue=0.001, outfmt=5)
            std_out, std_err = cline()

            # parse the output
            blast_records = NCBIXML.parse(StringIO(std_out))
            record = next(blast_records)

            # iterate over the test alignments
            # 1) verify the D3R criteria are satisfied
            # 2) create a test object & add it to the test_list
            for alignment in record.alignments:
                # check D3R criteria for each test alignment
                if not ParseAlignment.filter_blast_result(alignment, record):
                    continue
                else:
                    # extract pdb_id and a list of pdb_ids already assigned
                    pdb_id = ParseAlignment.get_alignment_pdb_id(alignment)

                    if not self.test_list:
                        self.instantiate_test(alignment,record)
                    else:
                        assigned = [x.get_pdb_id() for x in self.test_list]
                        # if the test object exists then append the sequence
                        if pdb_id in assigned:
                            chain_id = ParseAlignment.get_alignment_pdb_chain(alignment)
                            index = assigned.index(pdb_id)
                            self.test_list[index].set_sequence(pdb_id, chain_id)
                        # if it doesn't exist, instantiate the object and add it to the list
                        else:
                            self.instantiate_test(alignment,record)
            os.remove(fasta_file)
            #FastaIO.delete_fasta(fasta_file)

    def sort_by_resolution(self):
        """
        If a test list exists, and each test object has a resolution value, the list is sorted from highest resolution
        (lowest value) to lowest resolution (highest value)
        """
        self.test_list.sort(key=lambda x: x.get_resolution())

    def filter_by_experiment(self, method_type):
        """
        Removes test structures from test_list that weren't determined by method, method_type.
        :param method_type: (string) an experimental structure determination method used to solve wwpdb structures.
        Method types could include:
            'X-RAY DIFFRACTION', 'SOLUTION NMR', 'SOLID-STATE NMR', 'ELECTRON MICROSCOPY', 'ELECTRON CRYSTALLOGRAPHY',
            'FIBER DIFFRACTION', 'NEUTRON DIFFRACTION', or 'SOLUTION SCATTERING'.
        For the purpose of the CELPP cross-docking exercises, only structures determined by X-RAY DIFFRACTION are used.
        """
        del_test_indices = []
        assigned_test_structures = [x.get_pdb_id() for x in self.test_list]

        # identify test objects to delete
        for test in self.test_list:
            if test.get_expt_method() != method_type.upper():
                del_test_indices.append(assigned_test_structures.index(test.get_pdb_id()))

        # remove identified test objects
        del_test_indices.sort(reverse=True)
        for test_index in del_test_indices:
            del self.test_list[test_index]

    def instantiate_test(self, alignment, record):
        """
        Create a test instance, populate the attributes, and append it to the test list
        :param alignment:
        :param record
        """
        test = Test()

        pdb_id = ParseAlignment.get_alignment_pdb_id(alignment)
        chain_id = ParseAlignment.get_alignment_pdb_chain(alignment)

        # set the pdb_id
        test.set_pdb_id(pdb_id)

        # determine coverage & identity
        coverage, identity = ParseAlignment.coverage_and_identity(alignment, record)

        # set the test sequence
        test.set_sequence(pdb_id, chain_id)

        # set the coverage & identity values
        test.set_coverage(coverage)
        test.set_identity(identity)

        # get the experimental method
        test.set_expt_method()

        # parse the pdb and pull out the resolution
        test.set_resolution()

        # parse the pdb and pull ligand information
        test.set_ligands()

        self.test_list.append(test)

