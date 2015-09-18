__author__ = 'robswift'

import sys
from collections import defaultdict
from d3r.PreRelease.filtering_sets import do_not_call
from d3r.Blast.Target import Target

def add_ligands(non_polymer):
    """
    Read undesirable residues from the non_polymer_file, and return a list of target objects.
    :param non_polymer: absolute path to the pre-release non_polymer.tsv file
    :return: targets
    """
    ligands = read_n_label(non_polymer)
    ligands = remove_n_reduce(ligands)
    targets = instantiate_targets(ligands)


def read_n_label(non_polymer):
    """
    Reads ligands from new_release_nonpolymer.tsv. Each ligand is labeled 'do_not_call' or 'dock' depending on whether
    or not the ligand's resname is in the do_not_call filtering set. The ligand's resname, inchii string, and label are
    then added to the ligands dictionary, using the corresponding PDBID as a key. Specifically, a ligands dictionary
    entry has the following structure {'PDBID' :  [(resname, inchi, label), (resname, inchi, label), ... ]}, with each
    tuple representing a ligand associated with the PDB structure labeled by PDBID.
    :param non_polymer: the path to the new_release_nonpolymer.tsv file
    :return: ligand_dict, defaultdict(list)
    """
    ligands = defaultdict(list)
    handle = open(non_polymer, 'r')  # add a file checker class

    for line in handle.readlines()[1:]:
        words = line.split()
        if words:
            try:
                pdb_id = words[0].upper()
                resname = words[1]
                inchi = words[2]
                if resname and resname in do_not_call:
                    label = 'do_not_call'
                else:
                    label = 'dock'
                tup = (resname, inchi, label)
                ligands[pdb_id].append(tup)
            except:
                continue
    handle.close()
    return ligands


def remove_n_reduce(ligands):
    """
    Iterates over the contents for the ligands dictionary, and deletes structures if all of the ligands are labeled
    'do_not_call'
    :param ligands: {'PDBID' :  [(resname, inchi, label), (resname, inchi, label), ... ]}
    :return: ligands
    """
    for pdb in ligands.keys():
        lig_set = set([tup[0] for tup in ligands[pdb]])
        if lig_set.issubset(do_not_call):
            del ligands[pdb]
            # write to a log file
    return ligands

def instantiate_targets(ligands):
    """
    Creates target objects for each of the PDBIDs that are keys in the input ligands dictionary.
    :param ligands: {'PDBID' :  [(resname, inchi, label), (resname, inchi, label), ... ]}
    :return: targets, a list of target objects
    """
    targets = []
    for pdb in ligands.keys():
        target = Target()
        target.set_pdb_id(pdb)
        for ligand in ligands[pdb]:
            resname = ligand[0]
            inchi = ligand[1]
            label = ligand[2]
            target.set_ligand(format(resname), format(inchi), format(label))
        targets.append(target)
    return targets

def add_sequences(polymer, target_list):
    # target_list will be populated if the weekly release contains docking competent ligands. If the list is not
    # populated, then print an error and exit with error code 1
    if not target_list:
        print 'None of the ligands passed the \"do not call filter\"'
        sys.exit(1)
    else:
        file_handle = open(polymer,'r')
        target_ids = set([x.get_pdb_id() for x in target_list])

        for line in file_handle.readlines()[1:]:
            words = line.split()
            try:
                pdb_id = words[0].lower()
                chain_id = words[1]
                sequence = words[2]
                if pdb_id in target_ids:
                    index = [x.get_pdb_id() for x in target_list].index(pdb_id)
                    target_list[index].set_sequence(pdb_id, chain_id, sequence)
                else:
                    continue
            except:
                continue
        file_handle.close()
        return target_list
