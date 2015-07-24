__author__ = 'robswift'

from d3r.Blast.Test import Test

def get_alignment_pdb_id(alignment):
    """
    Returns a string of the four letter PDB ID found in alignment
    :param alignment:
    :return: 'pdb id'
    """
    pdb_chain_id = alignment.hit_def.encode('ascii').split()[0]
    pdb_id = pdb_chain_id.split('_')[0].lower()
    return pdb_id


def get_alignment_pdb_chain(alignment):
    """
    Returns a string of the chain id, e.g. 'A', 'B', etc.
    :param alignment:
    :return:
    """
    pdb_chain_id = alignment.hit_def.encode('ascii').split()[0]
    chain = pdb_chain_id.split('_')[1].upper()
    return chain


def filter_blast_result(alignment, record):
    """
    The binding pocket of the target structure should be covered by identical residues in the test sequence. We
    don't know where the target-sequence binding pocket is located. To increase the chances that the test sequence
    and target sequence have identical binding pockets, this function checks the following criteria:
    (i) X % of the target sequence must have a test-sequence alignment partner. Where X is given by the class
    variable, length threshold.
    (ii) Y % of the aligned sequence must be identical to the target sequence. Where Y is given by the class
    variable, identity_threshold.
    :param alignment: Bio.Blast.Record.Alignment object
    :param record: Bio.Blast.Record object
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


def blast_the_targets(target_list, pdb_db, fasta, out_dir, monomers=True):
    """
    Runs a blast search using the target sequence as a query. If monomers is True, then only monomers, or targets with
    one chain, will be used. If monomers is True, if no blast results that satisfy the monomer condition exist, then
    the target object is deleted from the target list.
    :param target_list: a list of target objects
    :param monomers: Boolean
    :return: target_list
    """
    Test.set_pdb_dict(fasta)

    if monomers:
        target_list = filter_by_chain(target_list, 1)
        if target_list:
            for target in target_list:
                target.run_blast(pdb_db, out_dir)
                target.filter_by_experiment('X-RAY DIFFRACTION')
                target.test_list = filter_by_chain(target.test_list, 1)
            target_list = clear_empty_targets(target_list)
    return target_list


def coverage_and_identity(alignment, record):
    aligned_test_seq_length = 0
    identical_test_seq_length = 0
    for hsp in alignment.hsps:
        aligned_test_seq_length += hsp.align_length
        identical_test_seq_length += hsp.identities
    coverage = float(aligned_test_seq_length) / record.query_length
    identity = float(identical_test_seq_length) / aligned_test_seq_length
    return coverage, identity


def filter_by_chain(molecule_list, chain_number):
    """
    Removes targets from target that contain a number of chains greater than chain_number
    :param molecule_list: may be either a list of target objects, or a list of test objects, e.g. target.test_list
    :param chain_number: an integer
    :return: target_list
    """
    if not molecule_list:
        return molecule_list
    else:
        tagged_for_removal = []
        assigned = [m.get_pdb_id() for m in molecule_list]

        for m in molecule_list:
            if len(m.sequences[m.pdb_id]) > chain_number:
                tagged_for_removal.append(assigned.index(m.pdb_id))

        tagged_for_removal.sort(reverse=True)
        for remove in tagged_for_removal:
            del molecule_list[remove]

    return molecule_list


def clear_empty_targets(target_list):
    """
    Removes target objects that are empty (they contain no test objects in test_list)
    :param target_list:
    :param target:
    :return:
    """
    tagged_for_removal = []
    assigned = [t.pdb_id for t in target_list]

    for target in target_list:
        if not target.test_list:
            tagged_for_removal.append(assigned.index(target.pdb_id))

    tagged_for_removal.sort(reverse=True)
    for remove in tagged_for_removal:
        del target_list[remove]

    return target_list


def remove_multiple_dockers(target_list):
    """
    Retains target structures that have only one dockable ligand
    :param target_list:
    :return:
    """
    tagged_for_removal = []
    assigned = [t.pdb_id for t in target_list]

    for target in target_list:
        if len(target.ligands) != 1:
            tagged_for_removal.append(assigned.index(target.pdb_id))

    tagged_for_removal.sort(reverse=True)
    for remove in tagged_for_removal:
        del target_list[remove]

    return target_list
