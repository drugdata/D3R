__author__ = 'robswift'

import csv
import os
import sys


def write_csv(target_list, out_dir):
    """

    :param target_list:
    :return:
    """
    if not os.access(out_dir, os.W_OK):
        print "%s is not writeable" % out_dir
        sys.exit(1)

    for target in target_list:
        target_pdb_id = target.get_pdb_id()
        file = os.path.join(out_dir, target_pdb_id + '.csv')

        out_handle = open(file, 'w')
        csv_writer = csv.writer(out_handle)

        # write target information
        csv_writer.writerow(['Target Information'])
        csv_writer.writerow(['PDBID', 'Ligand'])
        line = [target.get_pdb_id(), target.get_ligand_names()[0]]
        csv_writer.writerow(line)

        # write test information
        csv_writer.writerow([])
        csv_writer.writerow(['Test Information'])
        header_base = ['PDBID', 'Coverage', 'Identity', 'Resolution']
        max_no_ligands = max([len(test.get_ligand_names()) for test in target.test_list])
        if max_no_ligands == 0:
            header = header_base + ['Ligand']
        else:
            header = header_base + ['Ligand%s' % format(x+1) for x in range(max_no_ligands)]
        csv_writer.writerow(header)

        target.sort_by_resolution()
        for test in target.test_list:
            test_pdb_id = test.get_pdb_id()
            if test_pdb_id != target_pdb_id:
                coverage = test.get_coverage()
                identity = test.get_identity()
                resolution = test.get_resolution()
                line = [test_pdb_id, coverage, identity, resolution]
                ligands = [x for x in test.get_ligand_names()]
                if not ligands:
                    ligands.append('Apo')
                line = line + ligands
                csv_writer.writerow(line)

        out_handle.close()
