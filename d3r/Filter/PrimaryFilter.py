__author__ = 'robswift'

class PrimaryFilter(object):
    """
    Initiated with a list of target objects, this class can run a number of different filtering operations that will
    label the target objects 'tosser' if they fail to pass any one of the filtering criteria. After initializing and
    running the filtering operations, the method get_filtered_targets() is called, and the filtered target_list is
    returned. For example
        filter = filter(target_list)
        filter.filter_by_target_chains(5)
        filter.filter_by_experiment('X-RAY DIFFRACTION')
        filtered_list = filter.get_filtered_targets()
    """
    # default filtering values
    identity_threshold = 0.95
    coverage_threshold = 0.90
    chain_threshold = 1
    dockable_ligand_threshold = 1
    method = 'X-RAY DIFFRACTION'

    def __init__(self, targets):
        """
        The filter class is initiated with unique target structures (target objects) contained in the targets list. If
        the target list was filled by reading new_release_structure_sequence.tsv, and
        new_release_structure_nonpolymer.tsv, there is one target object for each of the unique wwPDB ids in
        new_release_sequence.tsv. And if their are corresponding ligands in new_release_structure_nonpolymer.tsv,
        ligand information is also included in the target object.
        :param targets: (list) a list of unique target objects
        """
        self.targets = targets

    def get_filtered_targets(self):
        """
        Returns the list of targets.
        :return: (list) targets, a list of filtered target objects
        """
        return self.targets

    def filter_by_inchi_error(self):
        """
        Remove target structures whose dockable ligand(s) produced an error on attempting to produce a rdkit.Chem.rdchem
        Mol object from the input InChI.
        """
        for target in self.targets:
            if target.inchi_error:
                target.set_reason(9)

    def filter_test_structures_by_identity(self, threshold = None):
        """
        Remove test structures whose sequences are insufficiently identical to the target sequence. This is measured by
        percent identity, where percent identity is the fraction of the BLAST alignment identical to the query. In the
        case of multi-chain BLAST hits, all of the chains must satisfy the threshold, or the hit is labeled a tosser.
        The default behavior removes those test structures whose sequences have a percent identity less than 95%. Each
        test structure with one or more sequences whose percent identity is less than the input threshold will be
        labeled 'tosser'. If all of the test structures of a given target are labeled 'tosser', then the target
        structure is also labeled 'tosser.'
        :param threshold: (float) percent identity threshold. Default = 0.95
        """
        if not threshold: threshold = PrimaryFilter.identity_threshold
        for target in self.targets:
            for test in target.test_list:
                if len([chain.identity for chain in test.sequences if chain.identity < threshold]) > 0:
                    test.set_reason(1)
            if len([test for test in target.test_list if test.triage]) == len(target.test_list):
                target.set_reason(8)


    def filter_test_structures_by_coverage(self, threshold = None):
        """
        Remove test structures whose sequences don't adequately cover the target sequence. This is measured by percent
        coverage, where percent coverage is the fraction of query sequence covered by the alignment. In the case of
        multi-chain BLAST hits, all of the chains must satisfy the threshold, or the hit is labeled a tosser. The
        default behavior removes those test structures whose sequences have a percent identity less than 90%. Each
        test structure with one or more sequences whose percent coverage is less than the input threshold will be
        labeled 'tosser'. If all of the test structures of a given target are labeled 'tosser', then the target is also
        labeled 'tosser.'
        :param threshold: (float) percent coverage threshold. Default = 0.90
        """
        if not threshold: threshold = PrimaryFilter.coverage_threshold
        for target in self.targets:
            for test in target.test_list:
                if len([chain for chain in test.sequences if chain.coverage < threshold]) > 0:
                    test.set_reason(2)
                    break
            if len([test for test in target.test_list if test.triage]) == len(target.test_list):
                target.set_reason(8)

    def filter_test_structures_with_too_many_chains(self, threshold = None):
        """
        Remove test structures with too many unique chains. The default behavior removes test structures with more than
        one unique chain. Each test object whose unique sequences number more than the input threshold will be labeled
        'tosser'. If all of the test objects of a given target are labeled 'tosser', then the target is also labeled
        'tosser.'
        :param threshold: (int) threshold value for the number of unique sequences in each test structure. Default = 1
        """
        if not threshold: threshold = PrimaryFilter.chain_threshold
        for target in self.targets:
            for test in target.test_list:
                if test.chain_count > threshold:
                    test.set_reason(3)
            if len([test for test in target.test_list if test.triage]) == len(target.test_list):
                target.set_reason(8)

    def filter_targets_with_too_many_dockable_ligands(self, threshold = None):
        """
        Remove target structures with too many dockable ligands. The default behavior removes target structures with
        more than one dockable ligand. Each target object whose number of dockable ligands is greater than the input
        threshold will be labeled'tosser'.
        :param threshold: (int) threshold value for the number of unique dockable ligands in each target. Default = 1.
        """
        if not threshold: threshold = PrimaryFilter.dockable_ligand_threshold
        for target in self.targets:
            if target.dock_count > threshold:
                target.set_reason(4)

    def filter_targets_with_only_apo_test_structures(self):
        """
        If a blast search of the wwPDB using the target sequence only returns apo test structures, then the target will
        be labeled 'tosser'.
        """
        for target in self.targets:
            if len([test for test in target.test_list if test.dock_count == 0]) == len(target.test_list):
                for test in target.test_list:
                    test.set_reason(10)
                target.set_reason(5)

    def filter_targets_without_dockable_ligands(self):
        """
        Remove target structures without dockable ligands. I.e. If an input target structure contains no dockable
        ligands, it is labeled 'tosser'.
        """
        for target in self.targets:
            if target.dock_count == 0:
                target.set_reason(6)

    def filter_test_structures_by_method(self, method=None):
        """
        Removes test structures from that weren't determined by the input method type.
        Method types can include:
            'X-RAY DIFFRACTION', 'SOLUTION NMR', 'SOLID-STATE NMR', 'ELECTRON MICROSCOPY', 'ELECTRON CRYSTALLOGRAPHY',
            'FIBER DIFFRACTION', 'NEUTRON DIFFRACTION', 'SOLUTION SCATTERING'.
        :param method_type: (string) an experimental structure determination method used to solve wwPDB structures.
        Default = 'X-RAY DIFFRACTION'
        """
        if not method: method = PrimaryFilter.method
        for target in self.targets:
            for test in target.test_list:
                if test.exp_method != method.upper():
                    test.set_reason(7)
            if len([test for test in target.test_list if test.triage]) == len(target.test_list):
                target.set_reason(8)

