__author__ = 'robswift'

import sys

from d3r.blast.mcss import MCSS
import logging

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    from rdkit import DataStructs
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit.Chem import Descriptors
except ImportError:
    logger.exception('Unable to import rdkit Ligand class will not work')


class Ligand(object):
    """

    """
    # timeout in seconds for rdkit to find maximum common substructure
    FINDMCS_TIMEOUT = 60
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
                    logger.exception('Caught exception trying to read ' + str(compinchi) + ' file')
                    continue
            handle.close()
        except IOError:
            logger.exception('Caught exception attempting to open or parse Components.inchi.ich file')
            sys.exit(1)
        if len(Ligand.inchi_component.keys()) < 20000:
            raise IOError('There is a problem with Components-inchi.ich.')

    def __init__(self, resname=None, inchi=None):
        self.resname = resname
        self.inchi = inchi
        #modified by sliu 08/04, here to add chain info for ligand 
        self.chain = None 
        self.label = None
        self.rd_mol = None              # a rd mol object that represents
        #modified by sliu 08/08 to add info for the ligand size
        self.size = None 
        self.heavy = None
        self.rot = None
        self.mcsss = []                 # A list of d3r.blast.Mcss objects
        #self.tanimoto = []              # A list of (tanimoto_score, reference_ligand_name)
        self.tanimoto = {}              # A dic of reference_ligand_name:tanimoto_score
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
            logger.error('No inchi string set')
            return False
        elif not inchi and self.inchi:
            inchi = self.inchi
        try:
            self.rd_mol = Chem.MolFromInchi(format(inchi), removeHs=False, sanitize=False, treatWarningAsError=True)
            return True
        except:
            logger.exception('Unable to create rdkt.Chem.rdmol Mol object')
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
            logger.exception('Unable to create rdkit.Chem.rdmol from resname')
            return False
    
    def set_heavy_size (self,):
        """
        Get number of heavy atoms from the rd_mol object associated with this ligand
        """
        if not self.rd_mol:
            print("Please set an rd mol object before calculating num of heavy atoms")
            logger.error('No rd mol set')
            return False
        try: 
            self.heavy = self.rd_mol.GetNumHeavyAtoms()
            return True
        except:
            logger.exception('Unable to get num of heavy atoms from rd mol object')
            return False

    def set_size (self,):
        """
        Get number of atoms from the rd_mol object associated with this ligand
        """
        if not self.rd_mol:
            print("Please set an rd mol object before calculating num of atoms")
            logger.error('No rd mol set')
            return False
        try:
            self.size = len(self.rd_mol.GetAtoms())
            return True
        except:
            logger.exception('Unable to get num of heavy atoms from rd mol object')
            return False

    def set_rot(self,):
        """
        set number of rotatable bond
        """
        if not self.rd_mol:
            print("Please set an rd mol object before calculating num of atoms")
            logger.error('No rd mol set')
            return False
        try:
            self.rot = Descriptors.NumRotatableBonds(self.rd_mol)
            return True
        except:
            logger.exception('Unable to set the num of rotatable bond from rd mol object')
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
            logger.debug('Trying to find MCS')

            res = rdFMCS.FindMCS([reference.rd_mol, self.rd_mol], timeout=Ligand.FINDMCS_TIMEOUT)
            mcss = Chem.MolFromSmarts(res.smartsString)
            return mcss
        except:
            logger.exception('Caught exception attempting to run rdkit FindMCS or MolFromSmarts')
            return None

    def calc_tanimoto(self, reference):
        """
        Determin the tanimoto similarity score based on the rd mol fingureprint
        """
        try:
            fps_ref = FingerprintMols.FingerprintMol(reference.rd_mol)
            fps_self = FingerprintMols.FingerprintMol(self.rd_mol)
            tanimoto = DataStructs.FingerprintSimilarity(fps_ref, fps_self)
            return tanimoto
        except:
            logger.exception('Caught exception attempting to run rdkit FingerprintMols.FingerprintMol or DataStructs.FingerprintSimilarity')
            return None
    def set_mcss(self, reference, mcss_mol, tanimoto_score):
        """
        Creates an Mcss object and appends it to the list, Ligand.mcsss. It's assumed that the input maximum common
        substructure was determined using the input reference ligand and the ligand object to which the method belongs.
        :param ref_name: d3r.blast.ligand object with the rd_mol attribute set
        :param mcss_mol: the maximum common substructure (rd mol object)
        """
        logger.debug('set_mcss()')
        mcss = MCSS(reference.resname, mcss_mol)
        #here mcss object is stored the mcss rd mol returned from mcss fucntion 
        mcss.set_size()
        #mcss.set_heavy_size()
        mcss.test = self.resname
        mcss.tanimoto = tanimoto_score
        self.mcsss.append(mcss)
    def set_tanimoto(self, reference, tanimoto_score):
        """
        store the tanimoto score and corresponding reference ligand resname
        """
        #self.tanimoto.append((tanimoto_score, reference.resname))
        self.tanimoto[reference.resname] = tanimoto_score
