__author__ = 'robswift'


import os
import logging
from d3r.filter.filter import QueryFilter
from d3r.filter.filter import HitFilter

logger = logging.getLogger(__name__)


class InputAnalysis(object):
    """

    """
    def __init__(self, targets):
        """

        :param targets:
        :return:
        """
        self.targets = targets
        self.number_of_structures = None            # total number of structures in targets
        self.number_of_complexes = None             # total number of structures with at least one bound ligand
        self.one_dock = None                        # total number of structures with only one dockable ligand
        self.at_least_one_dock = None               # total number of structures with at least one dockable ligand
        self.number_of_monomers = None              # total number of structures that are monomers
        self.number_of_monomers_one_dock = None     # total number of structures that are monomers with one dockable ->
                                                    # -> ligand.
        self.number_of_multimers = None             # total number of structures with more than one unique chain.
        self.number_of_multimers_one_dock = None    # total number of structures that are multimers with one dockable ->
                                                    # -> ligand.

    def set_number_of_structures(self):
        if self.number_of_structures:
            pass
        else:
            self.number_of_structures = len(self.targets)

    def set_number_of_complexes(self):
        if self.number_of_complexes:
            pass
        else:
            self.number_of_complexes = 0
            for target in self.targets:
                if int(target.ligand_count) > 0:
                    self.number_of_complexes += 1

    def set_one_dock(self):
        if self.one_dock:
            pass
        else:
            self.one_dock = 0
            for target in self.targets:
                if int(target.dock_count) == 1:
                    self.one_dock += 1

    def set_at_least_one_dock(self):
        if self.at_least_one_dock:
            pass
        else:
            self.at_least_one_dock = 0
            for target in self.targets:
                if int(target.dock_count) > 0:
                    self.at_least_one_dock += 1

    def set_number_of_monomers(self):
        if self.number_of_monomers:
            pass
        else:
            self.number_of_monomers = 0
            for target in self.targets:
                if int(target.sequence_count) == 1:
                    self.number_of_monomers += 1

    def set_number_of_monomers_one_dock(self):
        if self.number_of_monomers_one_dock:
            pass
        else:
            self.number_of_monomers_one_dock = 0
            for target in self.targets:
                if int(target.sequence_count) == 1 and int(target.dock_count) == 1:
                    self.number_of_monomers_one_dock += 1

    def set_number_of_multimers(self):
        if self.number_of_multimers:
            pass
        else:
            self.number_of_multimers = 0
            for target in self.targets:
                if int(target.sequence_count) > 1:
                    self.number_of_multimers += 1

    def set_number_of_multimers_one_dock(self):
        if self.number_of_multimers_one_dock:
            pass
        else:
            self.number_of_multimers_one_dock = 0
            for target in self.targets:
                if int(target.sequence_count) > 1 and int(target.dock_count) == 1:
                    self.number_of_multimers_one_dock += 1

    def print_to_standard_out(self):
        self.set_number_of_structures()
        self.set_number_of_complexes()
        self.set_at_least_one_dock()
        self.set_one_dock()
        self.set_number_of_monomers()
        self.set_number_of_monomers_one_dock()
        self.set_number_of_multimers()
        self.set_number_of_multimers_one_dock()
        print "The total number of structures is: %s" % self.number_of_structures
        print "The total number of complexes is: %s" % self.number_of_complexes
        print "The total number of complexes with at least one dockable ligand is: %s" % self.at_least_one_dock
        print "The total number of complexes with only one dockable ligand is: %s" % self.one_dock
        print "The total number of monomers is: %s" % self.number_of_monomers
        print "The total number of monomers with one dockable ligand is: %s" % self.number_of_monomers_one_dock
        print "The total number of multimers is: %s" % self.number_of_multimers
        print "The total number of multimers with one dockable ligand is: %s" % self.number_of_multimers_one_dock

    def print_to_file(self, out_dir):
        self.set_number_of_structures()
        self.set_number_of_complexes()
        self.set_at_least_one_dock()
        self.set_one_dock()
        self.set_number_of_monomers()
        self.set_number_of_monomers_one_dock()
        self.set_number_of_multimers()
        self.set_number_of_multimers_one_dock()
        logger.debug('Writing summary.txt file to ' + out_dir)
        handle = open(os.path.join(out_dir,'summary.txt'), 'w')
        out = ['INPUT SUMMARY']
        out.append('  entries:{no:>32}'.format(no = self.number_of_structures))
        out.append('  complexes:{no:>30}'.format(no = self.number_of_complexes))
        out.append('  dockable complexes:{no:>21}'.format(no = self.one_dock))
        out.append('  monomers:{no:>31}'.format(no = self.number_of_monomers))
        out.append('  dockable monomers:{no:>22}'.format(no = self.number_of_monomers_one_dock))
        out.append('  multimers:{no:>30}'.format(no = self.number_of_multimers))
        out.append('  dockable multimers:{no:>21}'.format(no = self.number_of_multimers_one_dock))
        out.append('')
        for l in out:
            logger.debug(l)
            handle.write("%s\n" % l)
        handle.close()

class OutputAnalysis(object):
    """

    """
    targets = {}                    # {'PDB ID' : (sequence_count, no. of blast hits, no. of candidates, no. elected,
                                    #              pdb_ids)

    def __init__(self, query = None):
        if query is not None:
            self.query = query

    def count_hits(self):
        return len(self.query.hits)

    def candidate_count(self):
        return len([hit for hit in self.query.hits if not hit.triage])

    def elected_count(self):
        return len([hit for hit in self.query.hits if hit.retain])

    def pdb_ids(self):
        return ','.join([hit.pdb_id for hit in self.query.hits if hit.retain])

    def set_target_dict(self):
        pdb_id = self.query.pdb_id
        seq_count = self.query.sequence_count
        hit_count = self.count_hits()
        candidate_count = self.candidate_count()
        elected_count = self.elected_count()
        pdb_ids = self.pdb_ids()
        OutputAnalysis.targets[pdb_id] = (seq_count, hit_count, candidate_count, elected_count, pdb_ids)

    def print_filter_criteria(self, out_dir):
        handle = open(os.path.join(out_dir, 'summary.txt'), 'a')
        out = ['FILTERING CRITERIA']
        out.append('  No. of query sequences           <=    {no}'.format(no = QueryFilter.sequence_threshold))
        out.append('  No. of dockable ligands           =    {no}'.format(no = QueryFilter.dockable_ligand_threshold))
        out.append('  Percent identity                 >=    {no:.2}'.format(no = HitFilter.identity_threshold))
        out.append('  Percent Coverage                 >=    {no:.2}'.format(no = HitFilter.coverage_threshold))
        out.append('  No. of hit sequences             <=    {no}'.format(no = HitFilter.sequence_threshold))
        out.append('  Structure determination method:        {method}'.format(method = HitFilter.method))
        out.append('')
        for l in out:
            logger.debug(l)
            handle.write("%s\n" %l)
        handle.close()

    def print_to_file(self, out_dir):
        handle = open(os.path.join(out_dir, 'summary.txt'), 'a')
        out = ['OUTPUT SUMMARY']
        out.append('  Targets found:{no:>26}'.format(no=len(OutputAnalysis.targets.keys())))
        for query in OutputAnalysis.targets.keys():
            sc, hc, cc, ec, pdb_ids = OutputAnalysis.targets[query]
            out.append('  Target: {id}|Sequences: {sc}|Hits: {hc}|'
                       'Candidates: {cc}|Elected:{ec}|PDBids: {pdbs}'.format(id=query, sc=sc, hc=hc, cc=cc, ec=ec,
                                                                             pdbs=pdb_ids))
        for l in out:
            logger.debug(l)
            handle.write("%s\n" %l)
        handle.close()