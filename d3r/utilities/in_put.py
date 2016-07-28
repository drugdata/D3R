__author__ = 'robswift'

import logging
from d3r.filter.filtering_sets import do_not_call
from d3r.blast.query import Query

logger = logging.getLogger(__name__)


def create_queries(polymer, non_polymer, ph):
    """
    Runs InPut.read_sequences and Input.read_ligands, which read in the sequence information contained in
    new_release_sequence.tsv and the ligand information contained in new_release_nonpolymer.tsv, respectively, and
    returns a list of d3r.blast.Target() objects.
    :param polymer: (string) the absolute path to new_release_sequence.tsv
    :param non_polymer: (string) the absolute path to new_release_nonpolymer.tsv
    :return: (d3r.blast.Target() object)
    """
    queries = read_sequences(polymer)
    if queries is not None:
        logger.debug('Found ' + str(len(queries)) + ' queries from ' +
                     polymer + ' file')

    read_ph(ph, queries)
    read_ligands(non_polymer, queries)
    return queries


def read_sequences(polymer):
    """
    Reads the information in a new_release_sequence.tsv file. A blast.Target() object is created for each unique wwPDB
    ID and appended to a list, which is returned.
    :param polymer: the absolute path to a new_release_sequence.tsv file
    :return: queries a list of query objects.
    """
    queries = []
    handle = open(polymer, 'r')
    for line in handle.readlines():
        if line.startswith('PDB_ID'):
            continue

        words = line.split()
        try:
            pdb_id = words[0].lower()
            chain_id = words[1]
            seq = words[2]
            queries = add_sequence(queries, seq, pdb_id, chain_id)
        except:
            logger.exception('Caught exception')
            continue
    handle.close()
    return queries


def add_sequence(queries, seq, pdb_id, chain_id):
    """
    Adds a FASTA sequence, with a specific chain and wwPDB ID, to the
    appropriate target object.
    :param queries: a list of target objects, which may be empty
    :param seq: a FASTA sequence. It will be converted to a Bio.SeqRecord
    :param pdb_id: The corresponding wwPDB id.
    :param chain_id: The corresponding chain id.
    """
    added = [query.pdb_id for query in queries if query]
    if pdb_id in added:
        queries[added.index(pdb_id)].set_sequence(chain_id, seq)
    else:
        query = Query()
        query.pdb_id = format(pdb_id.lower())
        query.set_sequence(chain_id, seq)
        queries.append(query)
    return queries


def read_ligands(non_polymer, queries):
    """
    Reads the information in a new_release_structure_nonpolymer.tsv file. The
    ligands are mapped to the appropriate
    target object by their wwPDB IDs. Morerover, each ligand is labeled as
    'do_not_call' or 'dock' depending on whether
    or not the resname of the ligand is found in the do_not_call set in the
    filter.filtering_sets module.
    :param non_polymer: absolute path to the pre-release non_polymer.tsv file
    :param queries, a list of target objects.
    """
    handle = open(non_polymer, 'r')
    for line in handle.readlines():
        if line.startswith('PDB_ID'):
            continue
        words = line.split()
        if words:
            try:
                pdb_id = words[0].lower()
                resname = words[1]
                inchi = words[2]
                ligand_label = label(resname)
                add_ligand(pdb_id, resname, inchi, ligand_label, queries)
            except:
                logger.exception('Caught exception')
                continue
    handle.close()


def add_ligand(pdb_id, resname, inchi, label, queries):
    """
    Adds the ligand, represented by its resname, inchi string, and label ( which can be 'dock' or 'do_not_call'), to
    the target with the corresponding wwPDB id.
    :param pdb_id: (string) wwPDBID
    :param resname: (string) the resname of the ligand.
    :param inchi: (string) the inchi string of the ligand
    :param queries: (list) a list of target objects
    :return:
    """
    added = [target.pdb_id for target in queries if target]
    queries[added.index(pdb_id)].set_ligand(resname, inchi, label)

def read_ph(ph, queries):
    """

    :param ph:
    :param queries:
    :return:
    """
    handle = open(ph, 'r')
    for line in handle.readlines():
        if line.startswith('PDB_ID'):
            continue
        words = line.split()
        if words:
            try:
                pdb_id = words[0].lower()
                exp_ph = words[1]
                add_ph(pdb_id, exp_ph, queries)
            except:
                logger.exception('Caught exception')
                continue
    handle.close()

def add_ph(pdb_id, exp_ph, queries):
    """

    :param pdb_id:
    :param exp_ph:
    :param queries:
    :return:
    """
    added = [query.pdb_id for query in queries if query]
    queries[added.index(pdb_id)].exp_ph = exp_ph

def label(resname):
    """
    Returns a label for the input resname. The label is either 'do_not_call' or 'dock', depending on whether the ligand
    is found in the do_not_call set in the filter.filtering_sets.
    :param resname:
    :return: ligand_label, a string
    """
    if resname and resname in do_not_call:
        ligand_label = 'do_not_call'
    else:
        ligand_label = 'dock'
    return ligand_label
