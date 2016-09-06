__author__ = 'robswift'

class MCSS(object):
    """
    Stores a maximum common substructure be
    """
    def __init__(self, reference, rd_mol):
        self.reference = reference      # The name of the reference molecule in the MCSS calculation, e.g. -->
                                        # --> Target.dock[0].resname (string)
        self.test = None                # The name of the test molecule in the MCSS calculation. e.g.
                                        # --> Test.dock[0].resname (string)
        self.rd_mol = rd_mol            # A rd mol object that represents the maximum common substructure between -->
                                        # --> the reference molecule and the molecule that contains the MCSS object.
        self.size = None                # The number of atoms in the MCSS (int)
        self.heavy = None
        self.tanimoto = None

    def set_size(self):
        """
        Determines the size, or atom count, of the MCSS. Note that the mol attribute must be set before the MCSS size
        can be set. If the mol object is not set, the method does nothing. Typical usage example:
                mcss = MCSS()
                mcss.reference = ligand.resname
                mcss.mol = Chem.MolFromSmarts(res.smartsString)
                mcss.set_size()
        """
        if not self.rd_mol:
            pass
        else:
            self.size = len(self.rd_mol.GetAtoms())
    #added by sliu 08/08 to add the info for heavy atoms
    def set_heavy_size(self):
        if not self.rd_mol:
            pass
        else:
            self.heavy = self.rd_mol.GetNumHeavyAtoms()
