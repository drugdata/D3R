__author__ = 'robswift'

class Filter(object):
    """
    Initiated with a list of target objects, this class can run a number of different filtering operations that will
    label the target objects 'tosser' if they fail to pass any one of the filtering criteria. After initializing and
    running the filtering operations, the method get_filtered_targets() is called, and the filtered target_list is
    returned. For example
        filter = Filter(target_list)
        filter.filter_by_target_chains(5)
        filter.filter_by_experiment('X-RAY DIFFRACTION')
        filtered_list = filter.get_filtered_targets()
    """
    def __init__(self, targets):
        """
        Filter the unique structures (target objects) in the targets list. If the target list was filled by reading
        new_release_structure_sequence.tsv and new_release_structure_nonpolymer.tsv, there is one target object for
        each of the unique wwPDB ids in new_release_sequence.tsv, and if their are corresponding ligands in
        new_release_structure_nonpolymer.tsv, ligand information is also included in the target object.
        :param targets: (list) a list of target objects
        """
        self.targets = targets

    def get_filtered_targets(self):
        """
        Returns the list of targets.
        :return:
        """
        return self.targets

    def filter_by_percent_identity(self, threshold):
        """
        Each test object whose percent identity is less than the input threshold, will be labeled 'tosser'.
        :param threshold: (float) percent identity threshold.
        """
        for target in self.targets:
            for test in target.test_list:
                if test.identity < threshold:
                    test.set_reason(1)

    def filter_by_percent_coverage(self, threshold):
        """
        Each test object whose percent coverage is less than the input threshold, will be labeled 'tosser'.
        :param threshold: (float) percent identity threshold.
        """
        for target in self.targets:
            for test in target.test_list:
                if test.coverage < threshold:
                    test.set_reason(2)
            if len([test for test in target.test_list if test.triage]) == len(target.test_list):
                target.set_reason(3)

    def filter_by_target_chain(self, threshold):
        """
        Each target object whose unique sequences number more than the input threshold, will be labeled 'tosser'.
        :param threshold: (int) threshold value for the number of unique sequences in each target structure.
        """
        for target in self.targets:
            if int(target.chain_count) > threshold:
                target.set_reason(3)

    def filter_by_test_chain(self, threshold):
        """
        Each test object whose unique sequences number more than the input threshold, will be labeled 'tosser'.
        :param threshold: (int) threshold value for the number of unique sequences in each test structure
        """
        for target in self.targets:
            for test in target.test_list:
                if test.chain_count > threshold:
                    test.set_reason(3)
            if len([test for test in target.test_list if test.triage]) == len(target.test_list):
                target.set_reason(3)

    def filter_by_number_dockable(self, threshold):
        """
        Each target object whose number of dockable ligands is greater than the input threshold will be labeled
        'tosser'.
        :param threshold: (int) threshold value for the number of unique dockable ligands in each target
        """
        for target in self.targets:
            if target.dock_count > threshold:
                target.set_reason(4)
