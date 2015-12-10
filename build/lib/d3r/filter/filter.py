__author__ = 'robswift'

from Bio import Alphabet


class BaseFilter(object):
    """

    """
    def __init__(self, query):
        """
        The filter class is initiated with unique target structures (target objects) contained in the targets list. If
        the target list was filled by reading new_release_structure_sequence.tsv, and
        new_release_structure_nonpolymer.tsv, there is one target object for each of the unique wwPDB ids in
        new_release_sequence.tsv. And if their are corresponding ligands in new_release_structure_nonpolymer.tsv,
        ligand information is also included in the target object.
        :param query: (list) a list of unique target objects
        """
        self.query = query


class QueryFilter(BaseFilter):
    """

    """
    sequence_threshold = 1
    dockable_ligand_threshold = 1

    def __init__(self, *args, **kwargs):
        super(QueryFilter, self).__init__(*args, **kwargs)

    def filter_by_inchi_error(self):
        """
        Remove target structures whose dockable ligand(s) produced an error on attempting to produce a rdkit.Chem.rdchem
        Mol object from the input InChI.
        """
        if self.query.inchi_error:
            self.query.set_reason(9)

    def filter_apo(self):
        """
        Remove queries structures without dockable ligands. I.e. If an input queries structure contains no dockable
        ligands, it is labeled 'tosser'.
        """
        if self.query.dock_count == 0:
            self.query.set_reason(6)

    def filter_by_dockable_ligand_count(self, threshold = None):
        """
        Remove query structures with too many dockable ligands. The default behavior removes query structures with
        more than one dockable ligand. Each query object whose number of dockable ligands is greater than the input
        threshold will be labeled'tosser'.
        :param threshold: (int) threshold value for the number of unique dockable ligands in each query. Default = 1.
        """
        if threshold == None:
            threshold = QueryFilter.dockable_ligand_threshold
            if self.query.dock_count > threshold:
                self.query.set_reason(4)

    def filter_by_sequence_count(self, threshold = None):
        if threshold == None:
            threshold = QueryFilter.sequence_threshold
            if self.query.sequence_count > threshold:
                self.query.set_reason(3)

    def filter_by_sequence_type(self):
        for seq_rec in self.query.sequences:
            if not Alphabet._verify_alphabet(seq_rec.seq):
                self.query.set_reason(11)
                break

class HitFilter(BaseFilter):
    """
    Initiated with a list of target objects, this class can run a number of different filtering operations that will
    label the Hit objects 'tosser' if they fail to pass any one of the filtering criteria. After initializing and
    running the filtering operations, the method get_filtered_query() is called, and the filtered target_list is
    returned. For example
        filter = filter(target_list)
        filter.filter_by_target_chains(5)
        filter.filter_by_experiment('X-RAY DIFFRACTION')
        filtered_list = filter.get_filtered_query()
    """
    # default filtering values
    identity_threshold = 0.95
    coverage_threshold = 0.90
    sequence_threshold = 4
    method = 'x-ray diffraction'

    def __init__(self, *args, **kwargs):
        super(HitFilter, self).__init__(*args, **kwargs)

    def filter_by_identity(self, threshold = None):
        """
        Remove hit structures whose sequences are insufficiently identical to the query sequence. This is measured by
        percent identity, where percent identity is the fraction of the BLAST alignment identical to the query. In the
        case of multi-chain BLAST hits, all of the chains must satisfy the threshold, or the hit is labeled a tosser.
        The default behavior removes those hit structures whose sequences have a percent identity less than 95%. Each
        hit structure with one or more sequences whose percent identity is less than the input threshold will be
        labeled 'tosser'. If all of the hit structures of a given target are labeled 'tosser', then the target
        structure is also labeled 'tosser.'
        :param threshold: (float) percent identity threshold. Default = 0.95
        """
        if threshold == None:
            threshold = HitFilter.identity_threshold
            for hit in self.query.hits:
                query_alignments = []
                for chain in [chain for chain in hit.sequences if chain.blast_hit]:
                    chain.sort_by_identity()
                    query_alignments.append(chain.query_alignments[0])
                if len([query_alignment for query_alignment in query_alignments if query_alignment.identity
                       < threshold]) > 0:
                    hit.set_reason(1)
            if len([hit for hit in self.query.hits if hit.triage]) == len(self.query.hits):
                self.query.set_reason(8)

    def filter_by_coverage(self, threshold = None):
        """
        Remove hit structures whose sequences don't adequately cover the query sequence. This is measured by percent
        coverage, where percent coverage is the fraction of query sequence covered by the alignment. In the case of
        multi-chain BLAST hits, all of the chains must satisfy the threshold, or the hit is labeled a tosser. The
        default behavior removes those hit structures whose sequences have a percent identity less than 90%. Each
        hit structure with one or more sequences whose percent coverage is less than the input threshold will be
        labeled 'tosser'. If all of the hit structures of a given query are labeled 'tosser', then the query is also
        labeled 'tosser.'
        :param threshold: (float) percent coverage threshold. Default = 0.90
        """

        # If the highest coverage for one or more chains is less than the threshold, then triage
        if threshold == None:
            threshold = HitFilter.coverage_threshold
            for hit in self.query.hits:
                query_alignments = []
                for chain in [chain for chain in hit.sequences if chain.blast_hit]:
                    chain.sort_by_coverage()
                    query_alignments.append(chain.query_alignments[0])
                if len([query_alignment for query_alignment in query_alignments if query_alignment.coverage
                        < threshold]) > 0:
                    hit.set_reason(2)
            if len([hit for hit in self.query.hits if hit.triage]) == len(self.query.hits):
                self.query.set_reason(8)

    def filter_by_sequence_count(self, threshold = None):
        """
        Remove hit structures with too many unique sequences. The default behavior removes hit structures with more than
        found unique sequences. Each hit object whose unique sequences number more than the input threshold will be
        labeled 'tosser'. If all of the hit objects of a given query are labeled 'tosser', then the query is also
        labeled 'tosser.'
        :param threshold: (int) threshold value for the number of unique sequences in each hit structure. Default = 4
        """
        if threshold == None:
            threshold = HitFilter.sequence_threshold
            for hit in self.query.hits:
                if hit.sequence_count > threshold:
                    hit.set_reason(3)
            if len([hit for hit in self.query.hits if hit.triage]) == len(self.query.hits):
                self.query.set_reason(8)

    def filter_apo(self):
        """
        If a blast search of the wwPDB using the query sequence only returns apo hit structures, then the query will
        be labeled 'tosser'.
        """
        if len([hit for hit in self.query.hits if hit.dock_count == 0]) == len(self.query.hits):
            for hit in self.query.hits:
                hit.set_reason(10)
            self.query.set_reason(5)

    def filter_by_method(self, method=None):
        """
        Removes hit structures that weren't determined by the input method type.
        Method types can include:
            'x-ray diffraction', 'solution nmr', 'solid-state nmr', 'electron microscopy', 'electron crystallography',
            'fiber diffraction', 'neutron diffraction', 'solution scattering'.
        :param method_type: (string) an experimental structure determination method used to solve wwPDB structures.
        Default = 'x-ray diffraction'
        """
        if not method: method = HitFilter.method
        for hit in self.query.hits:
            if hit.exp_method != method.lower():
                hit.set_reason(7)
        if len([hit for hit in self.query.hits if hit.triage]) == len(self.query.hits):
            self.query.set_reason(8)


class CandidateFilter(BaseFilter):
    """
    Filters the queries and hits that haven't been discarded by the PrimaryFilter. Specifically, for each remaining
    query object, the BLAST hits that satisfy the following criteria are identified
        1. The hit whose ligand is most similar to the ligand of the corresponding query
        2. The hit whose ligand is least similar to the ligand of the corresponding query
        3. The hit whose structure has the highest resolution
        4. The apo hit with the highest resolution
    """
    def __init__(self, *args, **kwargs):
        super(CandidateFilter, self).__init__(*args, **kwargs)

    def filter_for_most_similar(self):
        # sort the hits list by decreasing MCSS size and increasing resolution
        # the most similar will be at the front of the list, the least similar will be at the end of the list.
        hits = [hit for hit in self.query.hits if not hit.triage and hit.largest_mcss]
        if hits:
            hits.sort(key=lambda hit: (int(hit.largest_mcss.size), float(hit.resolution)), reverse = True)
            hits[0].set_retain_reason(1) # picks off the largest mcss with the highest resolution crystal structure


    def filter_for_least_similar(self):
        # sort the hits by increasing MCSS size and increasing resolution
        hits = [hit for hit in self.query.hits if not hit.triage and hit.smallest_mcss]
        if hits:
            hits.sort(key=lambda hit: (int(hit.smallest_mcss.size), float(hit.resolution)))
            hits[0].set_retain_reason(2) # picks off the smallest mcss with the highest resolution crystal structure
            #smallest = hits[0].smallest_mcss.size
            #for hit in hits:
            #    if hit.smallest_mcss.size > smallest:
            #        break
            #    else:
            #        hit.set_retain_reason(2)

    def filter_holo(self):
        hits = [hit for hit in self.query.hits if hit.resolution and not hit.triage and hit.dock_count > 0]
        if hits:
            hits.sort(key=lambda hit: hit.resolution)
            lowest = hits[0].resolution
            for hit in hits:
                if hit.resolution > lowest:
                    break
                else:
                    hit.set_retain_reason(3)

    def filter_apo(self):
        hits = [hit for hit in self.query.hits if hit.resolution and not hit.triage and hit.dock_count == 0]
        if hits:
            hits.sort(key=lambda hit: hit.resolution)
            lowest = hits[0].resolution
            for hit in hits:
                if hit.resolution > lowest:
                    break
                else:
                    hit.set_retain_reason(4)