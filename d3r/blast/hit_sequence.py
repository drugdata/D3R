__author__ = 'robswift'

import logging

logger = logging.getLogger(__name__)


class HitSequence(object):
    """

    """
    def __init__(self, seq_record):
        """
        Holds a sequence from a blast hit and information about it from one or more query alignments. Hit sequence
        information is stored in the following attributes:
            - seq_record        - A protein sequence (Bio.SeqRecord object)
            - hit_pdb_id        - pdb_id of the blast hit (string)
            - hit_chain_id      - chain_id of the blast hit (int)
            - hit_sequence_id   - the sequence id of the blast hit (int)
            - blast_hit         - True if the chain/sequence was returned from blast, false otherwise (boolean)
        Information about the alignment with a particular query is stored in a QueryAlignment object, which are stored
        in the following:
            - query_alignments  - A list of QueryAlignment objects (list)
        """
        self.seq_record = seq_record
        self.hit_pdb_id = None
        self.hit_chain_id = None
        self.hit_sequence_id = None
        self.blast_hit = None
        self.query_alignments = []

    def set_query_alignment(self, record, alignment):
        """
        returns True if a QueryAlignment object was successfully create, otherwise returns False
        """
        logger.debug('In seq_query_alignment()')
        try:
            qa = QueryAlignment()
            pdb_id, chain_id = record.query.encode('ascii').split()[0].split('_')
            query_length = record.query_length
            qa.query_pdb_id = pdb_id
            qa.query_chain_id = chain_id
            qa.query_length = query_length
            qa.alignment = alignment
            qa.set_coverage_and_identity()
            self.query_alignments.append(qa)
            return True
        except:
            logger.exception('Caught exception')
            return False

    def sort_by_coverage(self):
        self.query_alignments.sort(key=lambda qa: qa.coverage, reverse = True)

    def sort_by_identity(self):
        self.query_alignments.sort(key=lambda qa: qa.identity, reverse = True)


class QueryAlignment(object):
    """

    """
    def __init__(self):
        """
        Stores information about the alignment between a blast hit and a particular Query.
            - query_pdb_id      - pdb_id of the query (string)
            - query_chain_id    - chain_id of the query (string)
            - identity          - Percent identity (float)
            - coverage          - Percent coverage (float)
            - alignment         - An BLAST alignment (Bio.Record Alignment object)
            - query_length      - The length of the query sequence used to generate the input alignment (int)
        """
        self.query_id = None
        self.query_pdb_id = None
        self.query_chain_id = None
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
            logger.error("The TestSequence instance attributes, " +
                         "Ligand.query_length and Ligand.alignment " +
                         "both must be set " +
                         "before calling the set_coverage_and_identity method.")
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
        logger.debug('In _calculate_coverage_and_identity()')
        aligned_test_seq_length = 0
        identical_test_seq_length = 0
        for hsp in alignment.hsps:
            aligned_test_seq_length += hsp.align_length
            identical_test_seq_length += hsp.identities
        coverage = float(aligned_test_seq_length) / query_length
        identity = float(identical_test_seq_length) / aligned_test_seq_length
        return coverage, identity