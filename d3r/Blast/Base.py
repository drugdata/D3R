__author__ = 'robswift'


class Base(object):
    def __init__(self):
        self.pdb_id = None                  # 4 letter PDB code
        self.sequences = []                 # [Bio.SeqRecord, Bio.SeqRecord]
        self.chain_count = 0                # an integer count of the number of PDB chains
        self.dock = []                      # A list of dockable ligand objects.
        self.do_not_call = []               # A list of non-dockable ligand objects.
        self.dock_count = 0                 # an integer count of the number of dockable ligands
        self.ligand_count = 0               # an integer count of the total no. ligands
        self.triage = None                  # a string, 'tosser' or None, whose value determines whether a ->
                                            # -> structure (target or test object) should be discarded ('tosser') or ->
                                            # retained for subsequent docking (None).
        self.reason = None                  # a string that stores the reason, or failed filtering criteria, leading ->
                                            # -> to a triage decision of 'tosser'. The value of reason is selected ->
                                            # -> from the reasons dict by the function self.set_reason.
        self.reasons = {
            1 : 'The percent identity of the BLAST hit is below the threshold',
            2 : 'The percent coverage of the BLAST hit is below the threshold',
            3 : 'The number of unique sequences (chains) is above the threshold',
            4 : 'The number of dockable ligands is above the threshold',
            5 : 'A BLAST search of the wwPDB only returned apo structures',
            6 : 'No ligands are bound to the structure associated with the input pre-release sequence',
            7 : 'The structure of the BLAST hit was determined by a unsuitable method',
            8 : 'The BLAST query returned no suitable hits',
            9 : 'Error reading the InChI of the dockable ligand corresponding to the input pre-release sequence',
        }

    def set_reason(self, selection):
        """
        Sets the value of 'reason', an attribute that explains why the 'triage' attribute was set to 'tosser'. In
        addition to setting a reason, this method also sets the value of 'triage' to 'tosser'. To set a reason, an
        integer argument between 1 and 6, whose value maps to a particular reason, is provided. Selection to reason
        mapping:
            1 -> The percent identity of the BLAST hit is below the threshold
            2 -> The percent coverage of the BLAST hit is below the threshold
            3 -> The number of unique sequences (chains) is above the threshold
            4 -> The number of dockable ligands is above the threshold
            5 -> There are no ligand-bound structures that meet the BLAST criteria available in the wwPDB
            6 -> No ligands are bound to the structure associated with the input pre-release sequence
            7 -> The structure of the BLAST hit was determined by a unsuitable method
            8 -> The BLAST query returned no suitable hits
            9 -> Error reading the InChI of the dockable ligand corresponding to the input pre-release sequence
        Note that if the selection is outside of the range 1 to 9, no reason is set. Once the value for reason has
        been initiated, it is not fixed, and can be reset.
        :param selection: (int)
        """
        if int(selection) > 9 or int(selection) < 1:
            pass
        else:
            self.triage = 'tosser'
            self.reason = self.reasons[int(selection)]