__author__ = 'robswift'
__project__ = 'blastnfilter'

import os
from d3r.PreRelease import ParsePreRelease
from d3r.Blast import ParseAlignment
import OutPut


def run(options):
    non_polymer = os.path.abspath(options.non_polymer)
    polymer = os.path.abspath(options.polymer)
    out_dir = os.path.abspath(options.out)
    blast_dir = os.path.abspath(options.blast_db)
    pdb_db = os.path.join(blast_dir, 'pdb_db')
    fasta = os.path.join(blast_dir, 'pdb_seqres.txt')

    target_list = ParsePreRelease.add_ligands(non_polymer)
    target_list = ParsePreRelease.add_sequences(polymer, target_list)
    #new = [x for x in target_list if x.get_pdb_id().lower() == '2n02']
    target_list = ParseAlignment.blast_the_targets(target_list, pdb_db, fasta, out_dir)
    target_list = ParseAlignment.remove_multiple_dockers(target_list)
    OutPut.write_csv(target_list, out_dir)
