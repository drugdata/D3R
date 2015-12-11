__author__ = 'robswift'

import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS
from mcss import MCSS

class Ligand(object):
    """

    """

    inchi_component = {}

    @staticmethod
    def set_inchi_component(compinchi):
        """
        Each ligand resname is mapped to a InChi string. This information is stored in a class-wide dictionary
        inchi_component with the following structure
        inchi_component = { 'resname' : 'inchi' },
        where 'resname' is the resname used by the PDB for that ligand and 'inchi' is the corresponding PDB InChi
        string.
        :param compinchi: path to the wwPDB file Components-inchi.ich
        """
        try:
            handle = open(compinchi, 'r')
            for line in handle.readlines():
                try:
                    words = line.split()
                    inchi = words[0]
                    resname = words[1]
                    Ligand.inchi_component[format(resname)] = format(inchi)
                except:
                    continue
            handle.close()
        except IOError:
            sys.exit(1)
        if len(Ligand.inchi_component.keys()) < 20000:
            raise IOError('There is a problem with Components-inchi.ich.')

    def __init__(self, resname=None, inchi=None):
        self.resname = resname
        self.inchi = inchi
        self.label = None
        self.rd_mol = None              # a rd mol object that represents
        self.mcsss = []                 # A list of d3r.blast.Mcss objects
        self.mcss_error = None          # set to True, False or None

    def set_rd_mol_from_inchi(self, inchi = None):
        """
        Creates a rdkit.Chem.rdmol Mol object from an input InChI and sets it to self.rd_mol. Hydrogen atoms are not
        removed on successful conversion. If the conversion is unsuccessful, False is returned. If the conversion is
        successful, True is returned.
        :param inchi: (string) An InChI string
        :return: Boolean
        """
        if not inchi and not self.inchi:
            print("Please set an inchi string")
            return False
        elif not inchi and self.inchi:
            inchi = self.inchi
        try:
            self.rd_mol = Chem.MolFromInchi(format(inchi), removeHs=False, sanitize=False, treatWarningAsError=True)
            return True
        except:
            return False

    def set_rd_mol_from_resname(self, resname):
        """
        Creates a rdkit.Chem.rdmol from an input resname. If the rdmol is successfully created, the Ligand instance
        attribute, inchi is also set, and True is returned. Otherwise False is returned.
        :param resname: a wwPDB ID (string)
        :return: Boolean
        """
        try:
            inchi = Ligand.inchi_component[resname]
            self.rd_mol = Chem.MolFromInchi(format(inchi), removeHs=False, sanitize=False, treatWarningAsError=True)
            self.inchi = inchi
            return True
        except:
            return False

    def mcss(self, reference):
        """
        Determines the maximum common substructure (MCSS) between the input, or reference ligand and itself. The MCSS is
        stored in an RDkit mol object, which is returned. If no MCSS was found, the returned mol object will be empty.
        If an error prevents the MCSS calculation from completing, None is returned.
        Before this method is called, the ligand must have a rd_mol attribute set. For example
            lig = Ligand()
            lig.set_rd_mol_from_resname(resname) # where res is a wwPDB ID string
            mcss_mol = lig.mcss(reference)
            lig.set_mcss(reference.resname, mcss_mol)
        :param reference: d3r.blast.Ligand object with the rd_mol attribute set
        :return: the MCSS (rd mol object) or None
        """
        try:
            res = rdFMCS.FindMCS([reference.rd_mol, self.rd_mol])
            mcss = Chem.MolFromSmarts(res.smartsString)
            return mcss
        except:
            return None

    def set_mcss(self, reference, mcss_mol):
        """
        Creates an Mcss object and appends it to the list, Ligand.mcsss. It's assumed that the input maximum common
        substructure was determined using the input reference ligand and the ligand object to which the method belongs.
        :param ref_name: d3r.blast.ligand object with the rd_mol attribute set
        :param mcss_mol: the maximum common substructure (rd mol object)
        """
        mcss = MCSS(reference.resname, mcss_mol)
        mcss.set_size()
        mcss.test = self.resname
        self.mcsss.append(mcss)