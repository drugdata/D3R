__author__ = 'robswift'
__project__ = 'blastnfilter'

import os
import in_put
from d3r.blast.hit import Hit
from d3r.blast.ligand import Ligand
from d3r.filter.filter import QueryFilter
from d3r.filter.filter import HitFilter
from d3r.filter.filter import CandidateFilter
import out_put

def split_input(options):
    """
    Split out and return command line options
    :param options:
    :return: non_polymer (string), polymer (string), out_dir (string), blast_dir (string), pdb_db (string), ->
    -> fasta (string)
    """
    non_polymer = os.path.abspath(options.non_polymer)
    polymer = os.path.abspath(options.polymer)
    ph = os.path.abspath(options.ph)
    out_dir = os.path.abspath(options.out)
    blast_dir = os.path.abspath(options.blast_db)
    pdb_path = os.path.abspath(options.pdb_path)
    pdb_db = os.path.join(blast_dir, 'pdb_db')
    fasta = os.path.join(blast_dir, 'pdb_seqres.txt')
    compinchi = os.path.abspath(options.compinchi)
    return non_polymer, polymer, ph, out_dir, blast_dir, pdb_db, pdb_path, fasta, compinchi


def blast_the_query(query, pdb_db, pdb_path, fasta, out_dir, compinchi):
    """
    Runs a blast search using the target sequence as a query
    :param query: an instance of blast.query.Query
    :param pdb_db: The absolute path to a BLASTP database
    :param pdb_path: The absolute path to a decompressed copy of the PDB
    :param fasta: A file with fasta sequences of each chain in the PDB
    :param out_dir: The path to the directory where results will be written
    :param compinchi: The absolute path to a file with InChI strings for each ligand in the PDB
    :return: query
    """
    if not Hit.pdb_dir:
        Hit.set_pdb_dir(pdb_path)
    if not Hit.pdb_dict:
        Hit.set_pdb_dict(fasta)
    if not Ligand.inchi_component:
        Ligand.set_inchi_component(compinchi)
    records = query.run_blast(pdb_db, out_dir)
    if records:
        query.set_hits(records)
        query.fill_sequences()
    return query

def calculate_mcss(query):
    """
    Calculates the maximum common substructure (MCSS) between each dockable-query ligand and each dockable-BLAST-hit
    ligand
    :param queries: a list of Query objects
    """
    if query.dock_count == 0:
        pass
    else:
        for hit in [hit for hit in query.hits if hit.dock_count > 0]:
            for hit_ligand in hit.dock:
                for query_ligand in query.dock:
                    mcss_mol = hit_ligand.mcss(query_ligand)
                    if mcss_mol:
                        hit_ligand.set_mcss(query_ligand, mcss_mol)
            hit.set_maxmin_mcss()

def query_filter(query):
    """
    Runs query filtering operations
    :param query:
    :return:
    """
    filter = QueryFilter(query)
    filter.filter_apo()
    filter.filter_by_sequence_count()
    filter.filter_by_dockable_ligand_count()
    filter.filter_by_inchi_error()
    filter.filter_by_sequence_type()

def hit_filter(query):
    """
    Runs blast hit filtering operations
    :param query:
    :return:
    """
    filter = HitFilter(query)
    filter.filter_by_coverage()
    filter.filter_by_identity()
    filter.filter_by_sequence_count()
    filter.filter_apo()
    filter.filter_by_method()

def candidate_filter(queries):
    """
    Runs candidate filtering operations
    :param queries:
    :return:
    """
    filter = CandidateFilter(queries)
    filter.filter_holo()
    filter.filter_apo()
    filter.filter_for_most_similar()
    filter.filter_for_least_similar()

def run(options):
    """
    Run the BlastNFilter components
    :param options:
    """
    non_polymer, polymer, ph, out_dir, blast_dir, pdb_db, pdb_path, fasta, compinchi = split_input(options)
    queries = in_put.create_queries(polymer, non_polymer, ph)
    out_put.input_analysis(out_dir, queries)
    out_analysis = out_put.OutController()
    out_analysis.print_filter_criteria(out_dir)
    while queries:
        query = queries.pop()
        query_filter(query)
        if not query.triage:
            query = blast_the_query(query, pdb_db, pdb_path, fasta, out_dir, compinchi)
            calculate_mcss(query)
            hit_filter(query)
            candidate_filter(query)
        out_analysis.set_query(query)
        out_put.writer(out_dir, query, True)
    out_analysis.print_to_file(out_dir)

