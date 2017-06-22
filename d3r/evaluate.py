#! /usr/bin/env python

import os
from Bio import PDB
import commands
import time
import sys

try:
    from openeye.oechem import *
except ImportError as e:
    sys.stderr.write('Unable to import openeye, evaluate.py will '
                     'not work ' + str(e) + '\n')
import glob
import re
import pickle


def get_distance (pos1, pos2):
    _dist_0 = float(pos1.split(",")[0])-float(pos2.split(",")[0])
    _dist_1 = float(pos1.split(",")[1])-float(pos2.split(",")[1])
    _dist_2 = float(pos1.split(",")[2])-float(pos2.split(",")[2])
    return (_dist_0**2 + _dist_1**2 + _dist_2**2)**0.5

def get_center(ligand_pdb):
    xyz_lines = open(ligand_pdb, "r").readlines()
    multi_ligand = False
    atom_list = []
    x = y = z = 0
    for xyz_line in xyz_lines:
        if "HETATM" in xyz_line:
            #logging.debug("Check the get center of this protein: %s for this ligand: %s"%(protein_file, ligname))                                                                                                               
            atom_name = xyz_line[12:16]
            if atom_name in atom_list:
                multi_ligand = True
            else:
                atom_list.append(atom_name)
                try:
                    x += float(xyz_line[30:38])
                    y+= float(xyz_line[38:46])
                    z+= float(xyz_line[46:54])                                  
                except:                                                         
                    logging.debug("Fatal error: Cannot find the XYZ coordinate for this ligand:%s"%ligand_pdb)
                    return False                                                
    lig_center = "%8.3f, %8.3f, %8.3f"%(x/len(atom_list), y/len(atom_list), z/len(atom_list))
    logging.debug("Ligand center for this case:%s is %s"%(ligand_pdb, lig_center))        
    return lig_center                                                       

def extract_ligand_name(blastnfilter_result_txt):
    #extract the ligand name from blastnfilter result
    result_txt = open(blastnfilter_result_txt, "r")
    result_lines = result_txt.readlines()
    for result_line in result_lines:
        if "ligand," in result_line:
            ligand_name = result_line.split(",")[1].split("\n")[0].strip()
    return ligand_name

def extract_LMCSS_ligand_name (blastnfilter_result_txt):
    result_txt = open(blastnfilter_result_txt, "r")
    result_lines = result_txt.readlines()
    for result_line in result_lines:
        if "LMCSS," in result_line:
            ligand_name = result_line.split(",")[2].strip()
    return ligand_name

def mcs(ref_mol, fit_mol):
    #do the mcs search and return OEMatch object.
    #ignore Hydrogen 
    OESuppressHydrogens(fit_mol)
    OESuppressHydrogens(ref_mol)
    #set atom and bond expression                                       
    atomexpr = OEExprOpts_AtomicNumber
    bondexpr = 0
    #do the mcs search, using defined atom, bond expression options 
    mcss = OEMCSSearch( ref_mol, atomexpr, bondexpr, True)
    mcss.SetMCSFunc(OEMCSMaxAtomsCompleteCycles(1.5) )
    #create a new match object to store mcs search info.                
    new_match_list = []
    new_match_dic  = {}
    i = 0
    j = 0
    for match1 in mcss.Match(ref_mol):                                      
        i += 1                                                              
        #write out match1 molecule                                      
        mol1 = OEGraphMol()
        OESubsetMol(mol1,match1, True)                                      
        ofs1 = oemolostream("match1_%s.pdb"%i)                              
        OEWriteMolecule(ofs1, mol1)                                         
        ofs1.close()                                                        
        for match2 in mcss.Match(fit_mol):                                  
            j += 1                                                          
            check_list = []                                                 
            #write out match2 molecule                                      
            new_match = OEMatch()                                           
            mol2 = OEGraphMol()                                             
            OESubsetMol(mol2,match2, True)                                  
            ofs2 = oemolostream("match2_%s.pdb"%j)
            OEWriteMolecule(ofs2, mol2)                                     
            ofs2.close()                                                    
            for mp1, mp2 in zip(match1.GetAtoms(), match2.GetAtoms()):      
                ref_name = mp1.target.GetName().strip()                     
                fit_name = mp2.target.GetName().strip()                     
                new_match.AddPair (mp1.target, mp2.target)                  
            #store the match info
            new_match_list.append(new_match)                                
            new_match_dic[new_match] = (["match1_%s_vs_match2_%s"%(i,j), check_list ])
    return new_match_dic

def rmsd_mcss(ref_struc, fit_struc):
    #This function use the openeye mcss calculation to get the atom mapping first and then calculate RMSD, if multiple atom mapping is avaiable, the lowest RMSM will be returned 
    reffs = oemolistream()
    reffs.open(ref_struc)
    fitfs = oemolistream()
    fitfs.open(fit_struc)
    refmol = OEGraphMol()
    OEReadMolecule(reffs, refmol)
    for fitmol in fitfs.GetOEGraphMols():
        #get all possible matching 
        ss = mcs(refmol, fitmol,)
        mcss_rmsd_list = []
        match_info = []
        for mcss in ss.keys():
            #calculate the RMSD based on atom mapping
            mcss_rmsd = OERMSD(refmol, fitmol, mcss)
            mcss_rmsd_list.append(mcss_rmsd)
            match_info.append(ss[mcss])
    try:
        minimum_mcss = min(mcss_rmsd_list)
        return minimum_mcss
    except:
        return False

def wait_and_check (filename, timestep = 5, how_many_times = 100):                                                                       #add some relaxing time                                                                                                         
    count = 0 
    while (count < how_many_times):
        if not os.path.isfile(filename):
            time.sleep(timestep)
            count = count + 1
        else:
            return True
        return False

def align_protein (template_complex, input_complex, output_complex, timestep = 5, how_many_times = 100):
    #use schrodinger binding site alignment to get the aligned structure
    try:
        running_dir = os.getcwd()
        target_dir = os.path.dirname(os.path.abspath(output_complex))
        output_complex_filename = os.path.basename(output_complex)
        #change to the target dir since the align binding site only allow to run locally
        os.chdir(target_dir)
        commands.getoutput("$SCHRODINGER/utilities/align_binding_sites %s %s -o %s"%(template_complex, input_complex, output_complex_filename))
        os.chdir(running_dir)
        return wait_and_check(output_complex_filename, timestep = timestep, how_many_times = how_many_times)
    except Exception as ex:
        logging.exception("Could not align the protein %s onto template %s"%(input_complex, template_complex))

def whole_protein_align (template_complex, input_complex, output_complex, timestep = 5, how_many_times = 100):
    #use schrodinger whole protein alignment
    try:
        running_dir = os.getcwd()
        target_dir = os.path.dirname(os.path.abspath(output_complex))
        output_complex_filename = os.path.basename(output_complex)
        os.chdir(target_dir) 
        commands.getoutput("$SCHRODINGER/utilities/structalign %s %s"%(template_complex, input_complex))
        rotated_protein = "rot-" + input_complex
        if wait_and_check(rotated_protein, timestep = timestep, how_many_times = how_many_times):
            commands.getoutput("mv %s %s"%(rotated_protein, output_complex))  
            return True
        else:
            return False
    except Exception as ex:
        logging.exception("Could not align the protein %s onto template %s using the whole protein alignment"%(input_complex, template_complex))

def extract_ligand_from_complex (complex_pdb_file, ligand_pdb_file, ligand_info = "UNK-900"):
    #here the default ligand info is from schrodinger structure convert default ligand name. 
    complex_file = open(complex_pdb_file, "r")
    complex_lines = complex_file.readlines()
    complex_file.close()
    ligand_lines = []
    ligid = ligand_info.split("-")
    for complex_line in complex_lines:
        if (complex_line[17:20].strip()==ligid[0] and complex_line[22:26].strip()==ligid[1]):
            ligand_lines.append(complex_line)
    ligand_file = open(ligand_pdb_file, "w")
    ligand_file.writelines(ligand_lines)
    ligand_file.close()

def get_ligand_info_from_ligand_file (ligand_pdb_file):
    ligand_f = open(ligand_pdb_file, "r")
    ligand_lines = ligand_f.readlines()
    ligand_f.close()
    ligand_info = None
    for ligand_line in ligand_lines:
        ligand_title = ligand_line[17:20].strip()
        ligand_resnum = ligand_line[22:26].strip()
        new_ligand_info = ligand_title + "-" + ligand_resnum
        if not ligand_info:
            #get the first ligand info
            ligand_info = new_ligand_info
        else:
            if ligand_info == new_ligand_info:
                pass
            else:
                logging.info("Get multiple ligand info, need to check while will keep the first ligand info...")
                pass
    return ligand_info

def merge_two_pdb (receptor, ligand, complex_pdb):
    try:
        complex_lines = []
        f1 = open(receptor, "r")
        protein_lines = f1.readlines()
        f1.close()
        for p_line in protein_lines:
            if p_line [:6] not in ["CONECT", "ENDMDL" ] and p_line [:3] not in ["END"]:
                complex_lines.append(p_line)
        complex_lines.append("TER   \n")
        f2 = open(ligand, "r")
        ligand_lines = f2.readlines()
        f2.close()
        for l_line in ligand_lines:                               
            if l_line [:6] not in ["REMARK", "MODEL ", "CONECT", "ENDMDL"] and l_line not in ["END"]:
                complex_lines.append(l_line)
        f3 = open(complex_pdb, "w")
        f3.writelines(complex_lines)
        f3.close()
        return True
    except Exception as ex:
        logging.exception("The receptor %s and ligand %s cannot be merged into complex %s"%(receptor, ligand, complex_pdb))
        return False


def convert_ligand_format (input_ligand, output_ligand):
    try:
        commands.getoutput("$SCHRODINGER/utilities/structconvert %s %s"%(input_ligand, output_ligand))
        return True
    except Exception as ex:
        logging.exception("This ligand %s cannot be convertted to %s, need to check the format"%(input_ligand, output_ligand))
        return False

#def generate_ligand_and_receptor(complex_filename,ligand_filename, receptor_filename, ligand_info, rest_ligand_info):
def generate_ligand_and_receptor(complex_filename,ligand_filename, receptor_filename, ligand_info):
    complex_f = open(complex_filename, "r")
    complex_lines = complex_f.readlines()
    complex_f.close()
    ligand_lines = []
    receptor_lines = []
    for line in complex_lines:
        if 'ATOM' in line[:4] or 'HETATM' in line[:6]:
            ligand_name = line[17:20].strip()
            ligand_chain = line[20:22].strip()
            ligand_resnum = line[22:26].strip()
            logging.debug("The ligand info for the atom in the complex is %s-%s-%s"%(ligand_name, ligand_resnum, ligand_chain))
            if "%s-%s-%s"%(ligand_name, ligand_resnum, ligand_chain)  == ligand_info:
                logging.debug("find a ligand which fits the ligand info")
                ligand_lines.append(line)
                #receptor will have the target ligand 
                receptor_lines.append(line)
            else:
                #if not "%s-%s-%s"%(ligand_name, ligand_resnum, ligand_chain) in rest_ligand_info:
                if not 'HETATM' in line[:6]:
                    #receptor will append all atoms except the other ligand
                    receptor_lines.append(line)
                else:
                    pass
    receptor_f = open(receptor_filename, "w")
    receptor_f.writelines(receptor_lines)
    receptor_f.close()
    ligand_f = open(ligand_filename, "w")
    ligand_f.writelines(ligand_lines)
    ligand_f.close()

#create the crystal structure class
#input the ligand ID, get all structure from crystal file

class crystal (object):
    def __init__(self, crystal_file, ligand_name, pdbid):
        self._crystal = crystal_file
        self._ligand_name = ligand_name
        self._crystal_path = os.path.dirname(self._crystal)
        self._crystal_name = os.path.basename(self._crystal)
        parser = PDB.PDBParser()
        writer = PDB.PDBIO()
        self._biostruc = parser.get_structure (pdbid, crystal_file) 
    def get_ligand_info (self):
        #extract the target ligand information
        self._ligand_residue_info = []
        for _model in self._biostruc:
            for _chain in _model:
                for _residue in _chain:   
                    _residue_ligand_name = _residue.get_resname()
                    if _residue_ligand_name == self._ligand_name: 
                        _residue_info = "%s-%s-%s"%(_residue.get_resname(), _residue.get_id()[1], _chain.get_id())
                        if _residue_info not in self._ligand_residue_info:
                            self._ligand_residue_info.append(_residue_info)
                            logging.debug("Get the residue info :%s"%_residue_info)
                        else:
                            logging.debug("Got multiple ligand with the same residue info %s"%_residue_info)
                            pass
                    else:
                        pass
    def get_ligand_and_receptor_files(self):
        #get the ligand pdb file and receptor file per ligand 
        self._ligand = []
        self._receptor = []
        self._valid_ligand_info = []
        if len(self._ligand_residue_info) > 0:
            for _ligand_residue_info in self._ligand_residue_info:
                #get the ligand info for the rest of the ligand other than the target ligand and remove all of those ligand when generate the receptor
                _crystal_filename, _crystal_file_extension = os.path.splitext(self._crystal_name)
                _ligand_filename = _crystal_filename + "_" + _ligand_residue_info + "_" + "ligand" + _crystal_file_extension
                _receptor_filename = _crystal_filename + "_" + _ligand_residue_info + "_" + "receptor" + _crystal_file_extension  
                #store the ligand and receptor at the same path folder as the crystal path
                _ligand_file_full_path = os.path.join(self._crystal_path, _ligand_filename)
                _receptord_file_full_path = os.path.join(self._crystal_path, _receptor_filename)
                try:
                    generate_ligand_and_receptor (self._crystal, _ligand_file_full_path, _receptord_file_full_path, _ligand_residue_info )
                    logging.debug("Generating the ligand %s and receptor %s "%(_ligand_file_full_path, _receptord_file_full_path))
                    self._ligand.append(_ligand_file_full_path)
                    self._receptor.append(_receptord_file_full_path)
                    self._valid_ligand_info.append(_ligand_residue_info)
                except Exception as ex:
                    logging.exception("Failed to generate the ligand and receptor")
                    pass
        else:
            logging.info("There is no valid ligand in the crystal file, pass the ligand and receptor extraction step") 
            pass


class LMCSS (object):
    #initial LMCSS complex released to the participant, used to check the distance of LMCSS ligand with the crystal ligand
     
    def __init__(self, LMCSS_complex, ligand_info, crystal_obj):
        self._LMCSS = LMCSS_complex 
        self._LMCSS_path = os.path.dirname(self._LMCSS)
        self._LMCSS_name = os.path.basename(self._LMCSS)
        self._LMCSS_basename = os.path.splitext(self._LMCSS_name)[0]
         
        self._ligand = ligand_info
        self._crystal = crystal_obj
    def align_LMCSS_onto_crystal(self, time_check_frequence = 5, check_point_number = 100,):
        self._all_aligned_complex = {}
        self._all_aligned_ligand = {}
        if self._LMCSS:
            for _crystal_receptor_index, _crystal_receptor in enumerate(self._crystal._receptor):
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_receptor_index]
                logging.debug("Try to align the docked ligand onto this crystal receptor:%s"%(_crystal_receptor))
                self._aligned_complex_name = self._LMCSS_basename + "_aligned" + "_" + _crystal_ligand_info + ".pdb"
                self._aligned_ligand_name = self._LMCSS_basename + "_aligned_ligand" + "_" + _crystal_ligand_info + ".pdb"
                #need to name the aligned complex
                self._aligned_complex = os.path.join(self._LMCSS_path, self._aligned_complex_name)
                self._aligned_ligand = os.path.join(self._LMCSS_path, self._aligned_ligand_name)
                #do the alignment
                logging.debug("Try to align the original LMCSS %s onto the crystal receptor %s"%(self._LMCSS, _crystal_receptor))
                #try the binding site alignment first
                if not align_protein(_crystal_receptor, self._LMCSS, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                    #go to whole protein alignment
                    total_relaxing_time = time_check_frequence*check_point_number
                    logging.debug("The binding site alignment from %s onto %s didn't finish in %s second... Need to break"%(self._LMCSS, _crystal_receptor, total_relaxing_time))
                    if not whole_protein_align(_crystal_receptor, self._LMCSS, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                        total_relaxing_time = time_check_frequence*check_point_number
                        logging.info("The whole protein alignment from %s onto %s didn't finish in %s second... Need to break"%(self._LMCSS, _crystal_receptor, total_relaxing_time))
                    else:
                        self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                else:
                    self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                logging.debug("Successfully align %s onto %s and get the aligned structure %s"%(self._LMCSS, _crystal_receptor, self._aligned_complex))
                #extract the ligand from aligned complex
                try:
                    extract_ligand_from_complex(self._aligned_complex, self._aligned_ligand, ligand_info = self._ligand)
                    logging.debug("Successfully extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))
                    self._all_aligned_ligand[_crystal_ligand_info] = self._aligned_ligand
                except Exception as ex:
                    logging.exception("Cannot extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))
    def calculate_distance(self):
        #calculate the distance between LMCSS ligand and all possible crystal ligand 
        self._dis = {}
        self._min_dis = None 
        if self._all_aligned_ligand:
            for _crystal_ligand_index, _crystal_ligand in enumerate(self._crystal._ligand):
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_ligand_index]
                if _crystal_ligand_info in self._all_aligned_ligand:
                    _docked_ligand_aligned = self._all_aligned_ligand[_crystal_ligand_info]
                    try:
                        _crystal_ligand_center = get_center(_crystal_ligand)
                        _docked_ligand_aligned_center = get_center(_docked_ligand_aligned)
                        _distance = get_distance(_crystal_ligand_center, _docked_ligand_aligned_center)
                        self._dis[_crystal_ligand_info] = _distance
                        logging.debug("The distance for %s vs %s is %s"%(_crystal_ligand, _docked_ligand_aligned, _distance))
                    except Exception as ex:
                        logging.exception("The distance calculation for %s vs %s failed"%(_crystal_ligand, _docked_ligand_aligned))
                if self._dis:
                    self._min_dis = sorted(self._dis.values())[0]

class docked (object):
    #need to name the title as LMCSS, SMCSS etc
    #1, merge ligand and receptor
    #2, align the complex onto the crystal receptor
    #3, calculate the RMSD for each aligned ligand 
    def __init__(self, docked_ligand_filename, docked_receptor_filename, docked_category, crystal_obj):
        #load the ligand file
        self._docked_type = docked_category
        self._docked_ligand_mol = docked_ligand_filename
        self._docked_ligand_mol_path = os.path.dirname(self._docked_ligand_mol)
        self._docked_ligand_mol_name = os.path.basename(self._docked_ligand_mol)
        self._docked_ligand_mol_basename, self._docked_ligand_mol_extension = os.path.splitext(self._docked_ligand_mol_name)
        #load the receptor file
        self._docked_receptor_pdb = docked_receptor_filename
        self._docked_receptor_pdb_path = os.path.dirname(self._docked_receptor_pdb)
        self._docked_receptor_pdb_name = os.path.basename(self._docked_receptor_pdb)
        self._docked_receptor_pdb_basename, self._docked_receptor_pdb_extension = os.path.splitext(self._docked_receptor_pdb_name)
        #check 1, if the path of ligand and the path of receptor are the same
        
        self._docked_ligand_pdb = None
        self._crystal = crystal_obj
        if self._docked_ligand_mol_extension != ".mol":
            logging.debug ("This docked ligand %s is not in mol file format"%self._docked_ligand_mol)
        elif self._docked_receptor_pdb_extension != ".pdb":
            logging.debug ("This docked receptor %s is not in pdb file format"%self._docked_receptor_pdb)
        else:
            #convert to pdb file
            self._docked_ligand_pdb_name = self._docked_ligand_mol_basename + "_ligand.pdb"
            #make sure the pdb format file has the same path with the mol format file
            self._docked_ligand_pdb = os.path.join(self._docked_ligand_mol_path, self._docked_ligand_pdb_name) 
            #convert mol to pdb
            if not convert_ligand_format(self._docked_ligand_mol, self._docked_ligand_pdb):
                self._docked_ligand_pdb = None
            
    def create_complex (self):
        self._docked_complex = None
        if self._docked_ligand_pdb:
            #merge the receptor and the ligand together
            self._docked_complex_name = self._docked_ligand_pdb_name + "_complex.pdb"
            self._docked_complex = os.path.join(self._docked_ligand_mol_path, self._docked_complex_name)
            if not merge_two_pdb(self._docked_receptor_pdb, self._docked_ligand_pdb, self._docked_complex):
                self._docked_complex = None
    def align_complex_onto_crystal (self, time_check_frequence = 5, check_point_number = 100, sc_ligand_default = "UNK-900"):
        self._all_aligned_complex = {} 
        self._all_aligned_ligand = {} 
        if self._docked_complex:
            for _crystal_receptor_index, _crystal_receptor in enumerate(self._crystal._receptor):
                #_crystal_ligand = self._crystal._crystal_ligand[_crystal_receptor_index]
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_receptor_index]
                logging.debug("Try to align the docked ligand onto this crystal receptor:%s"%(_crystal_receptor))
                self._aligned_complex_name = self._docked_ligand_mol_basename + "_complex_aligned" + "_" + _crystal_ligand_info + ".pdb"
                self._aligned_ligand_name = self._docked_ligand_mol_basename + "_complex_aligned_ligand" + "_" + _crystal_ligand_info + ".pdb"
                self._aligned_complex = os.path.join(self._docked_ligand_mol_path, self._aligned_complex_name)
                self._aligned_ligand = os.path.join(self._docked_ligand_mol_path, self._aligned_ligand_name)
                #do the alignment
                logging.debug("Try to align the docked complex %s onto the crystal receptor %s"%(self._docked_complex, _crystal_receptor))
                #try the binding site alignment first
                if not align_protein(_crystal_receptor, self._docked_complex, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                    #go to whole protein alignment
                    total_relaxing_time = time_check_frequence*check_point_number
                    logging.debug("The binding site alignment from %s onto %s didn't finish in %s second... Need to break"%(self._docked_complex, _crystal_receptor, total_relaxing_time))
                    if not whole_protein_align(_crystal_receptor, self._docked_complex, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                        total_relaxing_time = time_check_frequence*check_point_number
                        logging.info("The whole protein alignment from %s onto %s didn't finish in %s second... Need to break"%(self._docked_complex, _crystal_receptor, total_relaxing_time))
                    else:
                        self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                else:
                    self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                logging.debug("Successfully align %s onto %s and get the aligned structure %s"%(self._docked_complex, _crystal_receptor, self._aligned_complex))
                #extract the ligand from aligned complex
                try:
                    extract_ligand_from_complex(self._aligned_complex, self._aligned_ligand, ligand_info = sc_ligand_default)
                    logging.debug("Successfully extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))
                    self._all_aligned_ligand[_crystal_ligand_info] = self._aligned_ligand
                except Exception as ex:
                    logging.exception("Cannot extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))
    def calculate_rmsd_and_distance (self):
        # 
        self._rmsd_dis = {}
        self._min_rmsd_dis = None
        if self._all_aligned_ligand:
            for _crystal_ligand_index, _crystal_ligand in enumerate(self._crystal._ligand):
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_ligand_index]
                if _crystal_ligand_info in self._all_aligned_ligand:
                    _docked_ligand_aligned = self._all_aligned_ligand[_crystal_ligand_info]
                    try:
                        #_rmsd, _distance = rmsd_mcss (_crystal_ligand, _docked_ligand_aligned)
                        _rmsd = rmsd_mcss (_crystal_ligand, _docked_ligand_aligned)
                        _crystal_ligand_center = get_center(_crystal_ligand)
                        _docked_ligand_aligned_center = get_center(_docked_ligand_aligned) 
                        _distance = get_distance(_crystal_ligand_center, _docked_ligand_aligned_center)
                        _all_mapping_files = glob.glob("match*.pdb")
                        for _mapping_file in _all_mapping_files:
                            os.remove(_mapping_file)
                        self._rmsd_dis[_crystal_ligand_info] = (_rmsd, _distance)
                        logging.debug("The rmsd and distance for %s vs %s is %s"%(_crystal_ligand, _docked_ligand_aligned, _rmsd))
                    except Exception as ex:
                        logging.exception("The rmsd calculation for %s vs %s failed"%(_crystal_ligand, _docked_ligand_aligned))
            #get the lowest RMSD and the corresponding distance
            if self._rmsd_dis:
                self._min_rmsd_dis = sorted(self._rmsd_dis.values())[0]

def get_all_docked_type(score_dic, docked_type = "LMCSS"):
    #from a score dic like score_dic[target_ID][docked_type]
    #get all docked_type value list for all target ID
    docked_type_value = []
    if score_dic:
        for target_ID in score_dic:
            if docked_type in score_dic[target_ID]:
                docked_type_value.append(score_dic[target_ID][docked_type])
    return docked_type_value


def clean_up_list_of_value(list_of_value):
    """Given list this method creates a
    new list removing any items that have
    value of None
    :returns: list of values in same order with None
              values removed
    """
    #clean up all None values
    new_list = []
    if list_of_value is None:
        return new_list

    for item in list_of_value:
        if item:
            new_list.append(item)
    return new_list


def calculate_average_min_max_median(full_list_of_value):
    """Given a list of values in `full_list_of_value`
    generate basic stats about values in list
    :returns: tuple of strings (average, min, max, median)
    """
    list_of_value = clean_up_list_of_value(full_list_of_value)
    if list_of_value is not None and len(list_of_value) > 0:
        count = len(list_of_value)
        average = "%-15.3f, " % (float(sum(list_of_value))/count)

        if count == 1:
            median = "%-15.3f, " % list_of_value[0]
        else:
            sorted_list = sorted(list_of_value)
            median = "%-15.3f, " % sorted_list[int(round(count/2))]

        minimum = "%-15.3f, " % min(list_of_value)
        maximum = "%-15.3f, " % max(list_of_value)
    else:
        average = minimum = maximum = median = "%-15s, " % " "

    return average, minimum, maximum, median
        
class data_container(object):

    def __init__(self, ):
        self._data = {}

    def register(self, target_ID, docked_type, value):
        #the docked type could be either "LMCSS" or "LMCSS_dis"
        # the value could be either rmsd or distance
        if value or value == 0.0:
            if target_ID not in self._data:
                self._data[target_ID] = {}
            if docked_type:
                if docked_type not in self._data[target_ID]:
                    self._data[target_ID][docked_type] = value

    def layout_pickle( self, pickle_filename = "RMSD.pickle"):
        self._pf = open(pickle_filename, "w")
        pickle.dump(self._data, self._pf)
        self._pf.close()

    def layout_plain (self,  plain_filename = "RMSD"):
        self._plain_filename = plain_filename
        self._combined_csv_filename = self._plain_filename + ".csv" 
        self._combined_txt_filename = self._plain_filename + ".txt" 
        self._LMCSS_list = get_all_docked_type(self._data, docked_type = "LMCSS") 
        self._SMCSS_list = get_all_docked_type(self._data, docked_type = "SMCSS") 
        self._hiResApo_list = get_all_docked_type(self._data, docked_type = "hiResApo") 
        self._hiResHolo_list = get_all_docked_type(self._data, docked_type = "hiResHolo") 
        self._hiTanimoto_list = get_all_docked_type(self._data, docked_type = "hiTanimoto") 

        self._whole_data_lines = ["%-20s%-15s, %-15s, %-15s, %-15s, %-15s, "
                                  "%-15s \n" % ("Target_PDBID,", "LMCSS", "SMCSS",
                                                "hiResApo", "hiResHolo", "hiTanimoto",
                                                "LMCSS_ori_distance")]

        self._whole_data_lines.append("\nSummary Statistics\n\n")

        #append the number of cases line
        self._whole_data_lines.append("%-20s%-15s, %-15s, %-15s, %-15s, "
                                      "%-15s, \n" % ("Number_of_cases,",
                                                     len(self._LMCSS_list),
                                                     len(self._SMCSS_list),
                                                     len(self._hiResApo_list),
                                                     len(self._hiResHolo_list),
                                                     len(self._hiTanimoto_list)))

        #append the average min max median, first
        (self._average_line, self._max_line,
         self._min_line, self._medium_line) = ("%-20s"%"Average,",
                                               "%-20s"%"Maximum,",
                                               "%-20s"%"Minimum,",
                                               "%-20s"%"Median,")

        for value_list in (self._LMCSS_list, self._SMCSS_list,
                           self._hiResApo_list, self._hiResHolo_list,
                           self._hiTanimoto_list):
            (average_score, min_score,
             max_score,
             medium_score) = calculate_average_min_max_median(value_list)

            self._average_line += average_score 
            self._max_line += max_score 
            self._min_line += min_score 
            self._medium_line += medium_score
        self._average_line += "\n" 
        self._max_line += "\n" 
        self._min_line += "\n" 
        self._medium_line += "\n" 
        self._whole_data_lines.append(self._average_line)
        self._whole_data_lines.append(self._max_line)
        self._whole_data_lines.append(self._min_line)
        self._whole_data_lines.append(self._medium_line)
        self._whole_data_lines.append("\nIndividual Results\n\n")
        #append the value for each PDBID 
        all_pro_type = ["LMCSS", "SMCSS", "hiResApo", "hiResHolo",
                        "hiTanimoto", "LMCSS_ori"]
        for target_ID in self._data:
            self._new_line_score = ""
            self._valid_line = False
            for docked_type in all_pro_type:
                if docked_type in self._data[target_ID]:
                    docked_type_dis = docked_type + "_dis"
                    if docked_type_dis in self._data[target_ID]:
                        self._new_line_score += "%-6.3f(%-6.3f) , " % (self._data[target_ID][docked_type],
                                                                       self._data[target_ID][docked_type_dis])
                    else:
                        if not docked_type == "LMCSS_ori":
                            self._new_line_score += "%-15.3f, " % (self._data[target_ID][docked_type])
                        else:
                            self._new_line_score += "(%-6.3f)       , " % (self._data[target_ID][docked_type])
                    self._valid_line = True
                else:
                    self._new_line_score += "%-15s, " % (" ")
            if self._valid_line:
                pdbid_full = "%s," % target_ID
                self._new_line = "%-20s" % (pdbid_full)
                self._new_line += self._new_line_score
                self._new_line += "\n"
                self._whole_data_lines.append(self._new_line)
        #write out csv file
        self._plain_csv_file = open(self._combined_csv_filename, "w")
        self._plain_csv_file.writelines(self._whole_data_lines)
        self._plain_csv_file.close()
        #write out txt file
        self._whole_data_lines_txt = [item.replace(",",
                                                   " ") for item in self._whole_data_lines]
        self._plain_txt_file = open(self._combined_txt_filename, "w")
        self._plain_txt_file.writelines(self._whole_data_lines_txt)
        self._plain_txt_file.close()
    
     
        
def calculate_rmsd(crystal_obj, docked_lig_mol, docked_receptor, docked_structure_type):
    docked_obj = docked(docked_lig_mol, docked_receptor, docked_structure_type, crystal_obj)
    docked_obj.create_complex()
    docked_obj.align_complex_onto_crystal(check_point_number = 10)
    docked_obj.calculate_rmsd_and_distance()
    return docked_obj

def create_crystal_obj (crystal_file, ligand_name, crystal_id):
    crystal_obj = crystal(crystal_file, ligand_name, crystal_id)
    crystal_obj.get_ligand_info()
    crystal_obj.get_ligand_and_receptor_files()
    return crystal_obj

def get_submitted_file_list (file_pattern_name = "docked.mol"):
    all_mol_files = glob.glob("*%s"%file_pattern_name)
    all_receptor_files = []
    if all_mol_files:
        for mol_file in all_mol_files:
            receptor_file = mol_file.replace('.mol','.pdb')
            all_receptor_files.append(receptor_file)
    return (all_mol_files, all_receptor_files)

def extract_crystal_file (crystal_pdbid, crystal_path):
    crystal_middle_name = crystal_pdbid[1:3]
    crystal_ent_file = "pdb" + crystal_pdbid + ".ent"
    crystal_pdbloc = os.path.join(crystal_path, crystal_middle_name, crystal_ent_file)
    if os.path.isfile(crystal_pdbloc):
        return crystal_pdbloc
    else:
        return False

def main_score(dock_dir, pdb_protein_path, evaluate_dir, blastnfilter_dir, challenge_data_path, update= True):
    #1 copy docked dir info into evaluate dir and create score for targe
    #2 create score dir and copy the ideal files into score dir
    #3 create submit obj and crystal obj in the crystal obj, and store the cyrstal obj directly 
    #4 after each calculation, append the crystal obj to the final result file in the upper level path
    #5 if for individual folder, the crystal obj is empty, also save it and skip it and go to the next case
    #1 switch to evaluation folder
    os.chdir(evaluate_dir)
    current_dir = os.getcwd()
    result_container = data_container()
    result_pickle = "RMSD.pickle"
    pickle_full_path = os.path.join(current_dir, result_pickle) 
    result_plain = "RMSD"
    plain_full_path = os.path.join(current_dir, result_plain)
    #get all target dir info from the docked dir
    target_dirs = list(os.walk(dock_dir))[0][1]
    for target_dir in target_dirs:
        #valid target must have at least one RMSD value return, else this target is invalid
        valid_target = False
        
        #copy the target dir to the evaluate dir
        commands.getoutput("cp -r %s/%s %s"%(dock_dir,target_dir, evaluate_dir))
        target_name = os.path.basename(target_dir)
        logging.info("=================We start to work on this target %s================="%target_name)
        #2 switch into the individual case
        os.chdir(target_name)
        all_mol_files, all_receptor_files = get_submitted_file_list()
        if not all_mol_files:
            logging.info("For this folder %s, there is no valid docked structure"%target_dir)
            result_container.register(target_dir, docked_type = None, value = None)
            ######need to register into the result pickle and txt, csv 
            os.chdir(current_dir)
            continue
        else:
            os.mkdir("score")
            #copy the crystal structure first
            crystal_ID = target_name
            #extract the crystal ligand name from blastnfilter result txt file
            blastnfilter_result_name = crystal_ID + ".txt"
            blastnfilter_txt = os.path.join(blastnfilter_dir, blastnfilter_result_name)
            ligand_name = extract_ligand_name(blastnfilter_txt) 
            LMCSS_ligand_name = extract_LMCSS_ligand_name(blastnfilter_txt)
            #get crystal structure location
            crystal_file = extract_crystal_file(crystal_ID, pdb_protein_path)
            if crystal_file:
                pdbid_local_path = os.getcwd()
                commands.getoutput("cp %s score/crystal.pdb"%(crystal_file))
                try:
                    os.chdir("score")
                    crystal_obj = create_crystal_obj("crystal.pdb", ligand_name, crystal_ID)
                    logging.info("\tSuccessfully create the crystal object")
                    #calculate the distance between the LMCSS ligand and crystal ligand
                    #get the file location
                    LMCSS_ori_file_path = os.path.join(challenge_data_path, crystal_ID)
                    LMCSS_complex_only_list = []
                    LMCSS_ligand_only_list = []
                    all_LMCSS_files = glob.glob("%s/LMCSS*.pdb"%LMCSS_ori_file_path)
                    for LMCSS_file in all_LMCSS_files:
                        if "-lig.pdb" not in os.path.basename(LMCSS_file):
                            LMCSS_complex_only_list.append(LMCSS_file)
                        else:
                            LMCSS_ligand_only_list.append(LMCSS_file)
                    if len(LMCSS_complex_only_list) == 1 and len(LMCSS_ligand_only_list) == 1:
                        LMCSS_complex = LMCSS_complex_only_list[0]
                        LMCSS_ligand = LMCSS_ligand_only_list[0]
                        ligand_info = get_ligand_info_from_ligand_file(LMCSS_ligand) 
                        LMCSS_local = os.path.basename(LMCSS_complex)
                        #copy locally
                        commands.getoutput("cp %s %s"%(LMCSS_complex, LMCSS_local))
                        #create the LMCSS obj
                        LMCSS_obj = LMCSS(LMCSS_local, ligand_info, crystal_obj)
                        #align the
                        LMCSS_obj.align_LMCSS_onto_crystal(check_point_number = 10)
                        #calculate the distance
                        LMCSS_obj.calculate_distance()
                        LMCSS_distance = LMCSS_obj._min_dis
                        result_container.register(target_dir, docked_type = "LMCSS_ori", value = LMCSS_distance)
                        logging.info("\tSuccessfully calculate the distance between original LMCSS ligand vs crstal ligand. Distance is %s"%LMCSS_distance)
                    else:
                        logging.info("There are %s original LMCSS complex files, which is abnormal and need to check"%len(LMCSS_complex_only_list))
                    os.chdir(pdbid_local_path)
                except Exception as ex:
                    result_container.register(target_dir, docked_type = None, value = None)
                    logging.exception("For this folder %s, could not create the crystal object"%target_dir)
                    os.chdir(current_dir)
                    continue
                #here get the crystal obj and loop all docked structure 
                for mol_index, potential_mol in enumerate (all_mol_files):              
                    potential_pdb = all_receptor_files[mol_index]                       
                    commands.getoutput("cp %s score"%potential_mol)                     
                    commands.getoutput("cp %s score"%potential_pdb)                     
                    #change path to local score folder
                    os.chdir("score")
                    local_potential_mol = os.path.basename(potential_mol)               
                    local_potential_pdb = os.path.basename(potential_pdb)
                    parsed_name = re.findall('([a-zA-Z0-9]+)-([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_docked.mol',local_potential_mol)
                    docked_structure_type = parsed_name[0][0]
                    try:
                        docked_obj = docked(local_potential_mol, local_potential_pdb, docked_structure_type, crystal_obj)
                        docked_obj.create_complex()                       
                        docked_obj.align_complex_onto_crystal(check_point_number = 10)
                        docked_obj.calculate_rmsd_and_distance()
                        (rmsd, dis) = docked_obj._min_rmsd_dis
                        docked_structure_type_dis = docked_structure_type + "_dis"
                        #rmsd = docked_obj._min_rmsd_dis[0]
                        valid_target = True
                        result_container.register(target_dir, docked_type = docked_structure_type, value = rmsd) 
                        result_container.register(target_dir, docked_type = docked_structure_type_dis, value = dis) 
                        logging.info("\tSuccessfully calculate the rmsd for this category %s, the rmsd is %s"%(docked_structure_type, rmsd))
                        #update the pickle and txt csv file if there is valid case found 
                        result_container.layout_pickle (pickle_filename = pickle_full_path)
                        result_container.layout_plain (plain_filename = plain_full_path)
                        os.chdir(pdbid_local_path)
                    except Exception as ex:
                        result_container.register(target_dir, docked_type = docked_structure_type, value = None) 
                        os.chdir(pdbid_local_path)
                        logging.exception("For this type of strcutre: %s, the rmsd could not be calculated"%docked_structure_type)
                        continue
                logging.info("++++++++++++++++Finished for this target %s+++++++++++++++++"%target_name)
                os.chdir(current_dir)
            else:
                result_container.register(target_dir, docked_type = None, value = None)
                logging.info("For this folder %s, there is no valid crystal structure"%target_dir)
                ######need to register into the result pickle and txt, csv 
                os.chdir(current_dir)
                continue
    #after loop all possible target, layout the pickle and plain files    
    result_container.layout_pickle (pickle_filename = pickle_full_path)
    result_container.layout_plain (plain_filename = plain_full_path)
            

if ("__main__" == __name__) :
    import logging
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-d", "--dockdir", metavar="PATH",
                  help="PATH where we could find the docking stage output")

    parser.add_argument("-b", "--blastnfilterdir", metavar="PATH",
                  help="PATH where we could find the blastnfilter stage output")

    parser.add_argument("-c", "--challengedir", metavar="PATH",
                  help="PATH where we could find the challengedata stage output")

    parser.add_argument("-o", "--outdir", metavar="PATH",
                  help="PATH where we will run the evaluate stage")

    parser.add_argument("-p", "--pdbdb", metavar="PATH",
                  help="PDB DATABANK which we will "
                       "get the crystal structure")
    logger = logging.getLogger()
    logging.basicConfig(format='%(asctime)s: %(message)s',
                    datefmt='%m/%d/%y %I:%M:%S', filename='final.log',
                    filemode='w', level=logging.INFO)
    opt = parser.parse_args()
    dockDir = opt.dockdir
    evaluateDir = opt.outdir
    pdbloc = opt.pdbdb
    blastnfilterDir = opt.blastnfilterdir
    challengedataDir = opt.challengedir
    # running under this dir
    running_dir = os.getcwd()
    pickle_result = main_score(dockDir, pdbloc, evaluateDir, blastnfilterDir, challengedataDir)
    # move the final log file to the result dir
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s" % (log_file_path, evaluateDir))
    
