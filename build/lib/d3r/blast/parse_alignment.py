__author__ = 'robswift'

from d3r.blast.hit import Hit

def filter_blast_result(alignment, record):
    """
    The binding pocket of the target structure should be covered by identical residues in the test sequence. We
    don't know where the target-sequence binding pocket is located. To increase the chances that the test sequence
    and target sequence have identical binding pockets, this function checks the following criteria:
    (i) X % of the target sequence must have a test-sequence alignment partner. Where X is given by the class
    variable, length threshold.
    (ii) Y % of the aligned sequence must be identical to the target sequence. Where Y is given by the class
    variable, identity_threshold.
    :param alignment: Bio.blast.Record.Alignment object
    :param record: Bio.blast.Record object
    :return: True if criteria are satisfied otherwise return False
    """
    length_threshold = 0.9
    identity_threshold = 0.95
    minimum_coverage = round(length_threshold * record.query_length)        # the smallest no. of target sequence
                                                                            # AA that align with the test sequence.
    minimum_identity = round(minimum_coverage * identity_threshold)         # the smallest no. of identical AA in
                                                                            # the alignment.
    # sum over the HSP (High Scoring Pairs of aligned test sequence residues)
    aligned_test_seq_length = 0
    identical_test_seq_length = 0
    for hsp in alignment.hsps:
        aligned_test_seq_length += hsp.align_length
        identical_test_seq_length += hsp.identities

    # check criteria
    if (minimum_coverage <= aligned_test_seq_length) & (minimum_identity <= identical_test_seq_length):
        return True
    else:
        return False