__author__ = 'robswift'

import sys
from collections import defaultdict
from filtering_sets import do_not_call
from d3r.Blast.Target import Target

def add_ligands(non_polymer):
    """
    Remove undesirable residues from the non_polymer_file, and return a list of target objects.
    :param non_polymer: absolute path to the pre-release non_polymer.tsv file
    :return: target_list
    """
    # stores a list of target objects
    target_list = []

    # load ligands in the default dict, ligand_dict[pdbid] = [(resname, inchi, label), (resname, inchi, label), ... ]
    # label is either 'do_not_call' or 'dock'
    ligand_dict = defaultdict(list)

    np_file_handle = open(non_polymer, 'r')  # add a file checker class

    for line in np_file_handle.readlines()[1:]:
        words = line.split()
        if words:
            try:
                pdb_id = words[0].upper()
                resname = words[1]
                inchi = words[2]
                if resname in do_not_call:
                    label = 'do_not_call'
                else:
                    label = 'dock'
                tup = (resname, inchi, label)
                ligand_dict[pdb_id].append(tup)
            except:
                continue
    np_file_handle.close()

    # iterate over ligand dictionary contents, & delete structures if all of the ligands are 'do_not_call'
    for pdb_id in ligand_dict.keys():
        lig_set = set([tup[0] for tup in ligand_dict[pdb_id]])
        if lig_set.issubset(do_not_call):
            continue
        else:
            target = Target()
            target.set_pdb_id(pdb_id)
            for ligand_tup in ligand_dict[pdb_id]:
                resname = ligand_tup[0]
                inchi = ligand_tup[1]
                label = ligand_tup[2]
                target.set_ligand(format(resname),format(inchi),format(label))
            target_list.append(target)
    return target_list


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
