__author__ = 'robswift'

class TestSequence(object):
    """

    """
    def __init__(self, seq_record):
        """
        Holds a sequence and information about it, particularly, information pertaining to CELPP BLAST alignments
        Attributes:
            - seq_record    - A protein sequence (Bio.SeqRecord object)
            - pdb_id        - A pdb_id (string)
            - chain_id      - A chain_id (string)
            - identity      - Percent identity (float)
            - coverage      - Percent coverage (float)
            - alignment     - An BLAST alignment (Bio.Record Alignment object)
            - query_length   - The length of the query sequence used to generate the input alignment (int)
        """
        self.seq_record = seq_record
        self.pdb_id = format(id)
        self.chain_id = None
        self.identity = None
        self.coverage = None
        self.alignment = None
        self.query_length = None

    def set_coverage_and_identity(self):
        """
        Determines the coverage and identity from the alignment and query length attributes and sets these attributes.
        Because these two attributes are prerequisites of the calculation, they must be set before the method call:
            d3r.blast.TestSequence(seq_record)
            seq.alignment = alignment
            seq.query_length
            seq.set_coverage_and_identity()
        If the calculation is successfully carried out, and the attributes successfully assigned, True is returned,
        otherwise, False is returned.
        :return: Boolean
        :return: Boolean
        """
        if not self.query_length or not self.alignment:
            print "The TestSequence instance attributes, Ligand.query_length and Ligand.alignment both must be set " \
                  "before calling the set_coverage_and_identity method."
            return False
        else:
            coverage, identity = self._calculate_coverage_and_identity(self.alignment, self.query_length)
            self.coverage = coverage
            self.identity = identity
            return True

    def _calculate_coverage_and_identity(self, alignment, query_length):
        """
        Calculates percent identity and percent coverage as fractions that are returned as floats. A input sequence to
        BLASTP is called a 'query' and a returned sequence is called a 'hit' sequence. Using these definitions, the
        percent identity is defined as the fraction of aligned hit residues identical to the query sequence, i.e.
        number_identical / hit_length. The percent coverage is defined as the ratio of the aligned hit length to the
        length of the query sequence, i.e. aligned_hit_length / query_length.
        :param alignment:
        :param record:
        :return: (float) coverage, (float) identity
        """
        aligned_test_seq_length = 0
        identical_test_seq_length = 0
        for hsp in alignment.hsps:
            aligned_test_seq_length += hsp.align_length
            identical_test_seq_length += hsp.identities
        coverage = float(aligned_test_seq_length) / query_length
        identity = float(identical_test_seq_length) / aligned_test_seq_length
        return coverage, identity