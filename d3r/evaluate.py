#! /usr/bin/env python

import os
from Bio import PDB
import commands
import time
import sys
import shutil
import logging
import glob
import re
import pickle
import json
from argparse import ArgumentParser

logger = logging.getLogger('d3r.evaluate')

try:
    from openeye.oechem import *
except ImportError as e:
    sys.stderr.write('Unable to import openeye, evaluate.py will '
                     'not work ' + str(e) + '\n')

# Name of log file that is written to by logging framework
FINAL_LOG = 'final.log'


def get_distance (pos1, pos2):
    """
    Given two X, Y, Z centers in format from output of get_center() function
    using standard three dimensional distance formula:

    square root (X^2 + Y^2 + Z^2)

    :param pos1: string of format 'X, Y, Z' Where X, Y, Z are numbers
    :param pos2: string of format 'X, Y, Z' Where X, Y, Z are numbers
    :return: distance as a float
    """
    _dist_0 = float(pos1.split(",")[0])-float(pos2.split(",")[0])
    _dist_1 = float(pos1.split(",")[1])-float(pos2.split(",")[1])
    _dist_2 = float(pos1.split(",")[2])-float(pos2.split(",")[2])
    return (_dist_0**2 + _dist_1**2 + _dist_2**2)**0.5


def get_center(ligand_pdb):
    """
    Parses a file that has atoms in format of:

    ATOM      1  N   SER A 183      35.409   3.264 119.011  1.00 44.03           N
    ATOM      2  CA  SER A 183      34.119   3.185 118.315  1.00 32.52           C
    HETATM 1900  C1  156 A 450      34.748  18.156  88.412  1.00 17.97           C
    HETATM 1901  C2  156 A 450      35.791  17.534  87.632  1.00 16.35           C

    and for lines starting with HETATM with unique column 2 (C1 row above)
    add column 6 (34.748 column above) to x, column 7 (18.156 column above) to y,
    and column 8 (88.412 column above)
    to z. These values would then be divided by number of HETATM lines found
    with unique column 2.

    :param ligand_pdb: path to file with data like above which corresponds to
                       atoms pertaining to receptor aka target and ligand
    :return: string of format 'X, Y, Z'
             where  X is sum of column 6 divided by number of atoms aka rows
             where  Y is sum of column 7 divided by number of atoms aka rows
             where  Z is sum of column 8 divided by number of atoms aka rows
             All numbers are set via format of %8.3f
    """
    xyz_lines = open(ligand_pdb, "r").readlines()
    multi_ligand = False
    atom_list = []
    x = y = z = 0
    for xyz_line in xyz_lines:
        if "HETATM" in xyz_line:
            #logger.debug("Check the get center of this protein: %s for this ligand: %s"%(protein_file, ligname))                                                                                                               
            atom_name = xyz_line[12:16]
            if atom_name in atom_list:
                # this variable multi_ligand is never used, not sure of its purpose
                multi_ligand = True
            else:
                atom_list.append(atom_name)
                try:
                    # Example line:
                    # HETATM 1900  C1  156 A 450      34.748  18.156  88.412  1.00 17.97

                    # 34.748 would be added to x if using example line above
                    x += float(xyz_line[30:38])

                    # 18.156 would be added to y if using example line above
                    y+= float(xyz_line[38:46])

                    # 88.412 would be added to z if using example line above
                    z+= float(xyz_line[46:54])                                  
                except:                                                         
                    logger.debug("Fatal error: Cannot find the XYZ coordinate for this ligand:%s"%ligand_pdb)
                    return False                                                
    lig_center = "%8.3f, %8.3f, %8.3f"%(x/len(atom_list), y/len(atom_list), z/len(atom_list))
    logger.debug("Ligand center for this case:%s is %s"%(ligand_pdb, lig_center))        
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


def wait_and_check (filename, timestep = 5, how_many_times = 100):
    """Using a while loop and sleep this function waits `how_many_times`
       for the file set by `filename` to appear on the file system
       sleeping for `timestep` seconds between retries.
       :param filename: This is the file the loop checks for existence
                        If set to None method will always return False
       :param timestep: Sleep in seconds between retries
       :param how_many_times: # of times to check, must be a value > 0
                              otherwise code will return False
    """
    if filename is None:
        logger.error('filename passed into wait_and_check is None')
        return False

    count = 0
    while count < how_many_times:
        if not os.path.isfile(filename):
            logger.debug('Sleeping ' + str(timestep) + ' waiting for file ' + filename + ' to appear')
            time.sleep(timestep)
            count = count + 1
        else:
            return True
    return False


def align_protein (template_complex, input_complex, output_complex, timestep = 5, how_many_times = 100):
    """use schrodinger binding site alignment to get the aligned structure

       specifically $SCHRODINGER/utilities/align_binding_sites
       which looks to perform pairwise superposition of multiple structures.
       First argument is reference structure.

       When align_binding_sites.py is run the following files are created, assuming
       these are first and 2nd arguments passed in:
       crystal-ligand-residueId-chainid_receptor.pdb
       LMCSS-<target id>-<candidate id>-<ligand id>.pd

       Created files:
       crystal_<ligand id>-<residue id>-<chain id>_receptor.csv
       crystal_<ligand id>-<residue id>-<chain id>_receptor-matrix.csv
       crystal_<ligand id>-<residue id>-<chain id>_receptor.log

       NOTE: if -LOCAL parameter is passed then these 2 files are also output
             crystal_<ligand id>-<residue id>-<chain id>_receptor-merged-input.maegz
             crystal_<ligand id>-<residue id>-<chain id>_receptor-align-initial.mae

       And the output_complex file is created.

       NOTE: It might be good to refactor and add -WAIT flag first to tell command to wait
             for job completion


       :param template_complex: path to file usually named crystal-ligand-residueId-chainid_receptor.pdb
       :param input_complex: path to file usually named LMCSS-<target id>-<candidate id>-<ligand id>.pdb
       :param output_complex: path to output file
    """
    try:
        running_dir = os.getcwd()
        target_dir = os.path.dirname(os.path.abspath(output_complex))
        output_complex_filename = os.path.basename(output_complex)
        #change to the target dir since the align binding site only allow to run locally
        os.chdir(target_dir)
        logger.debug("Running: PYTHONPATH= $SCHRODINGER/utilities/align_binding_sites " + template_complex +
                      " " + input_complex + " -o " + output_complex_filename)
        val = commands.getoutput("PYTHONPATH= $SCHRODINGER/utilities/align_binding_sites %s %s -o %s"%(template_complex, input_complex, output_complex_filename))
        logger.debug("output from align_binding_sites: " + str(val))
        os.chdir(running_dir)
        return wait_and_check(output_complex_filename, timestep = timestep, how_many_times = how_many_times)
    except Exception as ex:
        logger.exception("Could not align the protein %s onto template %s"%(input_complex, template_complex))


def whole_protein_align (template_complex, input_complex, output_complex, timestep = 5, how_many_times = 100):
    """use schrodinger whole structure alignment. This function calls structalign
       with template_complex as first argument and input_complex as second argument.
       The output always appears to be named:

       rot-<template_complex file name>
       rot-<input_complex file name>

       The above always seems to be written in the current working directory.
       This function takes the rot-<input_complex file name> and moves it to `output_complex`

    :param template_complex: path to file usually named crystal-ligand-residueId-chainid_receptor.pdb
    :param input_complex: path to file usually named LMCSS-<target id>-<candidate id>-<ligand id>.pdb
    :param output_complex: path to output file
    :param timestep:
    :param how_many_times:
    :return:
    """
    try:
        running_dir = os.getcwd()
        target_dir = os.path.dirname(os.path.abspath(output_complex))
        output_complex_filename = os.path.basename(output_complex)
        os.chdir(target_dir)
        val = commands.getoutput("PYTHONPATH= $SCHRODINGER/utilities/structalign %s %s"%(template_complex, input_complex))
        logger.debug('output from structalign: ' + str(val))
        rotated_protein = "rot-" + input_complex
        if wait_and_check(rotated_protein, timestep = timestep, how_many_times = how_many_times):
            commands.getoutput("mv %s %s"%(rotated_protein, output_complex))
            return True
        else:
            return False
    except Exception as ex:
        logger.exception("Could not align the protein %s onto template %s using the whole protein alignment"%(input_complex, template_complex))


def extract_ligand_from_complex (complex_pdb_file, ligand_pdb_file, ligand_info = "UNK-900"):
    """
      Reads pdb file `complex_pdb_file` and writes lines where [17:20] matches <ligand id>
      from `ligand_info` parameter and [22:26] matches <ligand residue id> from `ligand_info`
      parameter
    :param complex_pdb_file: Input pdb file to parse
    :param ligand_pdb_file:  Output pdb file to write lines matching <ligand id>
    :param ligand_info: String of format <ligand id>-<ligand residue id>
    :return:
    """
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
    """Given a LMCSS-<target>_<target>-<ligandid>-lig.pdb
       file which is a fragment of a pdb file with data
       that looks like this:
       HETATM 1899  C1  156 A 450      34.748  18.156  88.412  1.00 17.97           C
       HETATM 1900  C2  156 A 450      35.791  17.534  87.632  1.00 16.35           C
       HETATM 1901  C3  156 A 450      36.702  18.055  86.759  1.00 15.19           C

       The code then gets ligand id from [17:20] of line
       and ligand residue id from [22:26] which is 156 and
       450 respectively from the data above.
       The for loop logs messages if non matching
       ligand ids and residues are found, but otherwise
       continues.

       :returns: string <ligand id>-<ligand residue id>
    """
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
                # not sure why this is here, the log
                # message is not helpful and the first
                # encountered ligand info is what is kept
                logger.info("Get multiple ligand info, need to check while will keep the first ligand info...")
                pass
    return ligand_info


def merge_two_pdb (receptor, ligand, complex_pdb):
    """
    Merges two pdb files. Currently called from create_complex() method
    from docked() class.
    Function opens `receptor` file and grabs all lines of text that do NOT
    start with:

    CONECT
    ENDMDL
    END

    and writes them to output file `complex_pdb` in same order and appends this line
    with 3 spaces and a new line:

    TER

    Function then opens `ligand` file and grabs all lines of text that do NOT
    start with:

    REMARK
    MODEL
    CONECT
    ENDMDL

    and do NOT have the word END anywhere in the line. These lines are then appended
    to `complex_pdb` file

    :param receptor: A file path usually named <docked category>-<target id>-<candidate id>_docked.pdb
    :param ligand: A file path usually named <docked category>-<target id>-<candidate id>_docked_ligand.pdb
    :param complex_pdb: usually  <docked category>-<target id>-<candidate id>_docked_ligand_complex.pdb
    :return: True upon success or False if there is an Exception raised
    """
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
        logger.exception("The receptor %s and ligand %s cannot be merged into complex %s"%(receptor, ligand, complex_pdb))
        return False


def convert_ligand_format (input_ligand, output_ligand):
    """
    Runs schrodinger structconvert to convert `input_ligand` to `output_ligand`
    Invocation by docked() class uses this to convert .mol files to .pdb.
    Unfortunately this code uses commands.getoutput and does not actually check
    the exit code for failure.

    :param input_ligand: Input .mol file
    :param output_ligand: Output .pdb file
    :return: True upon success or False if there is an exception of any kind
    """
    try:
        commands.getoutput("PYTHONPATH= $SCHRODINGER/utilities/structconvert %s %s"%(input_ligand, output_ligand))
        return True
    except Exception as ex:
        logger.exception("This ligand %s cannot be convertted to %s, need to check the format"%(input_ligand, output_ligand))
        return False


def generate_ligand_and_receptor(complex_filename,ligand_filename, receptor_filename, ligand_info):
    """This function takes a pdb file `complex_filename` which would be the target ie 1fcz.pdb/ent
       and writes out 2 new pdb files, one containing just the ligand matching data in `ligand_info`

    :param complex_filename: input target pdb file
    :param ligand_filename: output file to write lines from `complex_filename` file starting with
                            ATOM or HETATM and with matching `ligand_info` data
    :param receptor_filename: output pdb file to write lines from `complex_filename` file starting with
                            ATOM in addition any lines written to ligand_filename will be written here
    :param ligand_info: string in format of <ligand name>-<ligand residue number>-<chain id>
                        in case of 1fcz the ligand was 156-450-A
    """
    complex_f = open(complex_filename, "r")
    complex_lines = complex_f.readlines()
    complex_f.close()
    ligand_lines = []
    receptor_lines = []
    for line in complex_lines:
        # Example line:
        # HETATM 1900  C1  156 A 450      34.748  18.156  88.412  1.00 17.97           C
        #
        # In above 17:20 is 156, 20:22 is A and 22:26 is 450
        #
        if 'ATOM' in line[:4] or 'HETATM' in line[:6]:
            ligand_name = line[17:20].strip()
            ligand_chain = line[20:22].strip()
            ligand_resnum = line[22:26].strip()
            logger.debug("The ligand info for the atom in the complex is %s-%s-%s"%(ligand_name, ligand_resnum, ligand_chain))
            if "%s-%s-%s"%(ligand_name, ligand_resnum, ligand_chain)  == ligand_info:
                logger.debug("find a ligand which fits the ligand info")
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


class crystal (object):
    """Represents a crystal
    """

    def __init__(self, crystal_file, ligand_name, pdbid):
        """

        :param crystal_file: file name only invocation seens uses crystal.pdb
        :param ligand_name: name of ligand ie 156. For standard, it was obtained
                            from
                            stage.3.blastnfilter/1fcz.txt from ligand, line.
        :param pdbid: target_name for gold standard its 1fcz
        """
        self._crystal = crystal_file
        self._ligand_name = ligand_name
        # in only invocation seen just crystal.pdb is passed in so
        # code assumes file is in working directory and path will
        # be set to ''
        #
        self._crystal_path = os.path.dirname(self._crystal)

        # since crystal.pdb is passed in this will be unchanged
        # and set to crystal.pdb
        #
        self._crystal_name = os.path.basename(self._crystal)
        parser = PDB.PDBParser()

        # writer is never used
        writer = PDB.PDBIO()
        self._biostruc = parser.get_structure (pdbid, crystal_file)

    def get_ligand_info (self):
        """
          This method builds an internal list of residues from
          the crystal file passed in via the constructor using
          Biopython's PDB.PDBParser(). The list only contains
          entries whose residue name matches ligand_name passed
          in via the constructor.

          The contents of the list is are strings of this format:
          <residue name>-<residue id>-<chain id>

          Duplicate entries are ignored except for a message to debug
          of logger
        """
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
                            logger.debug("Get the residue info :%s"%_residue_info)
                        else:
                            logger.debug("Got multiple ligand with the same residue info %s"%_residue_info)
                            pass
                    else:
                        pass

    def get_ligand_and_receptor_files(self):
        #get the ligand pdb file and receptor file per ligand 
        self._ligand = []

        # self._receptor is directly read by LMCSS.align_LMCSS_onto_crystal
        # and by docked.align_complex_onto_crystal
        #
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
                    # generate_ligand_and_receptor
                    # writes out crystal_ligand-residueId-chainid_ligand.pdb
                    #            crystal_ligand-residueId-chainid_receptor.pdb
                    # which contain just the atoms for the ligand (in _ligand file) and atoms
                    # for both target and ligand in the _receptor file)
                    generate_ligand_and_receptor (self._crystal, _ligand_file_full_path, _receptord_file_full_path, _ligand_residue_info )
                    logger.debug("Generating the ligand %s and receptor %s "%(_ligand_file_full_path, _receptord_file_full_path))
                    self._ligand.append(_ligand_file_full_path)
                    self._receptor.append(_receptord_file_full_path)
                    self._valid_ligand_info.append(_ligand_residue_info)
                except Exception as ex:
                    logger.exception("Failed to generate the ligand and receptor")
                    pass
        else:
            logger.info("There is no valid ligand in the crystal file, pass the ligand and receptor extraction step") 
            pass


class LMCSS (object):
    """initial LMCSS complex released to the participant, used to check
     the distance of LMCSS ligand with the crystal ligand
    """
     
    def __init__(self, LMCSS_complex, ligand_info, crystal_obj):
        """
        Constructor.
        :param LMCSS_complex: Path to LMCSS-<target id>-<candidate id>-<ligand id>.pdb which
                              contains pdb for both ligand and I think candidate
        :param ligand_info: Should be a string in format <ligand id>-<ligand residue id>
        :param crystal_obj: Instance of crystal class loaded with <target id>
        """
        self._LMCSS = LMCSS_complex 
        self._LMCSS_path = os.path.dirname(self._LMCSS)
        self._LMCSS_name = os.path.basename(self._LMCSS)

        # this variable is LMCSS_complex file name only with .pdb suffix stripped off
        self._LMCSS_basename = os.path.splitext(self._LMCSS_name)[0]
         
        self._ligand = ligand_info
        self._crystal = crystal_obj

    def align_LMCSS_onto_crystal(self, time_check_frequence = 5, check_point_number = 100,):
        self._all_aligned_complex = {}
        self._all_aligned_ligand = {}
        if self._LMCSS:

            # self._crystal._receptor contains a list of file paths to files
            # with name crystal-ligand-residueId-chainid_receptor.pdb
            # which contain atoms for target and ligand
            #
            for _crystal_receptor_index, _crystal_receptor in enumerate(self._crystal._receptor):

                # self._crystal._valid_ligand_info is a list filled with
                # strings of format <residue name>-<residue id>-<chain id>
                # and is filled by crystal.get_ligand_info() and
                # by crystal.get_ligand_and_receptor_files()
                # the enumeration requires the order to be consistent
                # with self._crystal._receptor
                #
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_receptor_index]
                logger.debug("Try to align the docked ligand onto this crystal receptor:%s"%(_crystal_receptor))


                # _aligned_complex_name is set to
                # LMCSS-<target id>-<candidate id>-<ligand id>_aligned_<residue name>-<residue id>-<chain id>.pdb
                #
                self._aligned_complex_name = self._LMCSS_basename + "_aligned" + "_" + _crystal_ligand_info + ".pdb"

                # _aligned_ligand_name is set to
                # LMCSS-<target id>-<candidate id>-<ligand id>_aligned_ligand_<residue name>-<residue id>-<chain id>.pdb
                #
                self._aligned_ligand_name = self._LMCSS_basename + "_aligned_ligand" + "_" + _crystal_ligand_info + ".pdb"

                #
                # prepend directory which from code looks to be outdir/<target>/score/
                #
                self._aligned_complex = os.path.join(self._LMCSS_path, self._aligned_complex_name)
                self._aligned_ligand = os.path.join(self._LMCSS_path, self._aligned_ligand_name)

                #do the alignment
                logger.debug("Try to align the original LMCSS %s onto the crystal receptor %s"%(self._LMCSS, _crystal_receptor))
                #try the binding site alignment first
                # below align_binding_sites is tried and then structalign
                #
                if not align_protein(_crystal_receptor, self._LMCSS, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                    #go to whole protein alignment
                    total_relaxing_time = time_check_frequence*check_point_number
                    logger.debug("The binding site alignment from %s onto %s didn't finish in %s second... Need to break"%(self._LMCSS, _crystal_receptor, total_relaxing_time))
                    if not whole_protein_align(_crystal_receptor, self._LMCSS, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                        total_relaxing_time = time_check_frequence*check_point_number
                        logger.info("The whole protein alignment from %s onto %s didn't finish in %s second... Need to break"%(self._LMCSS, _crystal_receptor, total_relaxing_time))
                    else:
                        # _crystal_ligand_info is a string of format
                        # <residue name>-<residue id>-<chain id>
                        # and its being used as a key to map to the aligned output
                        # LMCSS-<target id>-<candidate id>-<ligand id>_aligned_<residue name>-<residue id>-<chain id>.pdb
                        #
                        self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                else:
                    # _crystal_ligand_info is a string of format
                    # <residue name>-<residue id>-<chain id>
                    # and its being used as a key to map to the aligned output
                    # LMCSS-<target id>-<candidate id>-<ligand id>_aligned_<residue name>-<residue id>-<chain id>.pdb
                    #
                    self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                logger.debug("Successfully align %s onto %s and get the aligned structure %s"%(self._LMCSS, _crystal_receptor, self._aligned_complex))
                #extract the ligand from aligned complex
                try:
                    # extract ligand from output pdb and save to self._aligned_ligand file name
                    # LMCSS-<target id>-<candidate id>-<ligand id>_aligned_ligand_<residue name>-<residue id>-<chain id>.pdb
                    extract_ligand_from_complex(self._aligned_complex, self._aligned_ligand, ligand_info = self._ligand)
                    logger.debug("Successfully extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))
                    self._all_aligned_ligand[_crystal_ligand_info] = self._aligned_ligand
                except Exception as ex:
                    logger.exception("Cannot extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))

    def calculate_distance(self):
        #calculate the distance between LMCSS ligand and all possible crystal ligand
        # contains a dictionary of
        # <residue name>-<residue id>-<chain id> => <distance between centers>
        self._dis = {}

        # contains smallest distance found in self._dis
        # and is referenced in main_score() around line 1133
        #
        self._min_dis = None 
        if self._all_aligned_ligand:
            # self._crystal._ligand contains a list of file paths to files
            # with name crystal-ligand-residueId-chainid_ligand.pdb
            # which contain atoms for ligand
            #
            for _crystal_ligand_index, _crystal_ligand in enumerate(self._crystal._ligand):
                # self._crystal._valid_ligand_info is a list filled with
                # strings of format <residue name>-<residue id>-<chain id>
                # and is filled by crystal.get_ligand_info() and
                # by crystal.get_ligand_and_receptor_files()
                # the enumeration requires the order to be consistent
                # with self._crystal._receptor
                #
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_ligand_index]
                if _crystal_ligand_info in self._all_aligned_ligand:
                    _docked_ligand_aligned = self._all_aligned_ligand[_crystal_ligand_info]
                    try:
                        _crystal_ligand_center = get_center(_crystal_ligand)
                        _docked_ligand_aligned_center = get_center(_docked_ligand_aligned)
                        _distance = get_distance(_crystal_ligand_center, _docked_ligand_aligned_center)
                        # build a map of <residue name>-<residue id>-<chain id> => <distance between centers>
                        self._dis[_crystal_ligand_info] = _distance

                        logger.debug("The distance for %s vs %s is %s"%(_crystal_ligand, _docked_ligand_aligned, _distance))
                    except Exception as ex:
                        logger.exception("The distance calculation for %s vs %s failed"%(_crystal_ligand, _docked_ligand_aligned))

                # I'm guessing this can be taken out of this for loop and calculated just once
                if self._dis:
                    self._min_dis = sorted(self._dis.values())[0]


class docked (object):
    #need to name the title as LMCSS, SMCSS etc
    #1, merge ligand and receptor
    #2, align the complex onto the crystal receptor
    #3, calculate the RMSD for each aligned ligand 
    def __init__(self, docked_ligand_filename, docked_receptor_filename, docked_category, crystal_obj):
        """
        Constructor
        :param docked_ligand_filename: filename of user submission of mol file
                                       usually in format like so;
                                       <docked category>-<target id>-<candidate id>_docked.mol
        :param docked_receptor_filename: filename of user submission in pdb file
                                       usually in format like so;
                                       <docked category>-<target id>-<candidate id>_docked.pdb
        :param docked_category: string usually set to one of following;
                                LMCSS, SMCSS, hiResApo, hiResHolo
        :param crystal_obj: instance of crystal class properly loaded with target
        """

        # sets _docked_type to LMCSS, SMCSS, hiResApo, or hiResHolo
        self._docked_type = docked_category

        # set to <docked category>-<target id>-<candidate id>_docked.mol
        self._docked_ligand_mol = docked_ligand_filename

        # directory where mol file resides, but this will be empty since
        # code calling this strips off the path before calling this
        # constructor
        self._docked_ligand_mol_path = os.path.dirname(self._docked_ligand_mol)

        # just name of file, but will be same cause code calling this already
        # strips off path
        self._docked_ligand_mol_name = os.path.basename(self._docked_ligand_mol)

        # strip off extension and put base name in self._docked_ligand_mol_basename and put
        # extension in self._docked_ligand_mol_extension
        self._docked_ligand_mol_basename, self._docked_ligand_mol_extension = os.path.splitext(self._docked_ligand_mol_name)

        # set to <docked category>-<target id>-<candidate id>_docked.pdb
        self._docked_receptor_pdb = docked_receptor_filename

        # set to empty string cause main_score strips off path already
        self._docked_receptor_pdb_path = os.path.dirname(self._docked_receptor_pdb)

        # set to same as self._docked_receptor_pdb
        self._docked_receptor_pdb_name = os.path.basename(self._docked_receptor_pdb)

        # strip of extension and put base name into self._docked_receptor_pdb_basename
        # and put extension in self._docked_receptor_pdb_extension
        self._docked_receptor_pdb_basename, self._docked_receptor_pdb_extension = os.path.splitext(self._docked_receptor_pdb_name)
        #check 1, if the path of ligand and the path of receptor are the same
        
        self._docked_ligand_pdb = None

        # set to crystal class instance
        self._crystal = crystal_obj

        # code merely logs if extensions are incorrect
        # with no real checks
        if self._docked_ligand_mol_extension != ".mol":
            logger.debug ("This docked ligand %s is not in mol file format"%self._docked_ligand_mol)
        elif self._docked_receptor_pdb_extension != ".pdb":
            logger.debug ("This docked receptor %s is not in pdb file format"%self._docked_receptor_pdb)
        else:
            #
            # convert to pdb file using schroding convert_ligand_format
            # by creating new file name by appending _ligand.pdb to basename of mol
            # file which is <docked category>-<target id>-<candidate id>_docked
            #
            self._docked_ligand_pdb_name = self._docked_ligand_mol_basename + "_ligand.pdb"

            # make sure the pdb format file has the same path with the mol format file
            self._docked_ligand_pdb = os.path.join(self._docked_ligand_mol_path, self._docked_ligand_pdb_name)

            # convert mol to pdb and only set self._docked_ligand_pdb to None if there is a
            # failure.
            #
            if not convert_ligand_format(self._docked_ligand_mol, self._docked_ligand_pdb):
                self._docked_ligand_pdb = None
            
    def create_complex (self):
        """
        Merges creates <docked category>-<target id>-<candidate id>_docked_ligand_complex.pdb
        by merging <docked category>-<target id>-<candidate id>_docked_ligand.pdb with
        docked category>-<target id>-<candidate id>_docked.pdb
        """
        self._docked_complex = None
        if self._docked_ligand_pdb:
            #merge the receptor and the ligand together
            self._docked_complex_name = self._docked_ligand_pdb_name + "_complex.pdb"
            self._docked_complex = os.path.join(self._docked_ligand_mol_path, self._docked_complex_name)
            if not merge_two_pdb(self._docked_receptor_pdb, self._docked_ligand_pdb, self._docked_complex):
                self._docked_complex = None

    def align_complex_onto_crystal(self, time_check_frequence=5,
                                   check_point_number = 100,
                                   sc_ligand_default = "UNK-900"):
        self._all_aligned_complex = {} 
        self._all_aligned_ligand = {} 
        if self._docked_complex:
            # self._crystal._receptor contains a list of file paths to files
            # with name crystal_ligand-residueId-chainid_receptor.pdb
            # which contain atoms for target and ligand
            #
            for _crystal_receptor_index, _crystal_receptor in enumerate(self._crystal._receptor):
                # self._crystal._valid_ligand_info is a list filled with
                # strings of format <residue name>-<residue id>-<chain id>
                # and is filled by crystal.get_ligand_info() and
                # by crystal.get_ligand_and_receptor_files()
                # the enumeration requires the order to be consistent
                # with self._crystal._receptor and self._crystal._ligand
                #
                _crystal_ligand_info = self._crystal._valid_ligand_info[_crystal_receptor_index]
                logger.debug("Try to align the docked ligand onto this crystal receptor:%s"%(_crystal_receptor))

                # set to <docked category>-<target id>-<candidate id>_docked_complex_aligned_<residue name>-<residue id>-<chain id>.pdb
                self._aligned_complex_name = self._docked_ligand_mol_basename + "_complex_aligned" + "_" + _crystal_ligand_info + ".pdb"

                # set to <docked category>-<target id>-<candidate id>_docked_complex_aligned_ligand_<residue name>-<residue id>-<chain id>.pdb
                self._aligned_ligand_name = self._docked_ligand_mol_basename + "_complex_aligned_ligand" + "_" + _crystal_ligand_info + ".pdb"

                # prepend path to self._aligned_complex_name
                self._aligned_complex = os.path.join(self._docked_ligand_mol_path, self._aligned_complex_name)

                # prepend path to self._aligned_ligand_name
                self._aligned_ligand = os.path.join(self._docked_ligand_mol_path, self._aligned_ligand_name)

                #do the alignment
                logger.debug("Try to align the docked complex %s onto the crystal receptor %s"%(self._docked_complex, _crystal_receptor))

                # try the binding site alignment first with
                # _crystal_receptor set to crystal_ligand-residueId-chainid_receptor.pdb
                # self._docked_complex set to <docked category>-<target id>-<candidate id>_docked_ligand_complex.pdb
                # self._aligned_complex set to <docked category>-<target id>-<candidate id>_docked_complex_aligned_<residue name>-<residue id>-<chain id>.pdb
                # if align_protein() fails try whole_protein_align()
                #
                if not align_protein(_crystal_receptor, self._docked_complex, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                    #go to whole protein alignment
                    total_relaxing_time = time_check_frequence*check_point_number
                    logger.debug("The binding site alignment from %s onto %s didn't finish in %s second... Need to break"%(self._docked_complex, _crystal_receptor, total_relaxing_time))
                    if not whole_protein_align(_crystal_receptor, self._docked_complex, self._aligned_complex, timestep = time_check_frequence, how_many_times = check_point_number ):
                        total_relaxing_time = time_check_frequence*check_point_number
                        logger.info("The whole protein alignment from %s onto %s didn't finish in %s second... Need to break"%(self._docked_complex, _crystal_receptor, total_relaxing_time))
                    else:
                        self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex
                else:
                    self._all_aligned_complex[_crystal_ligand_info] = self._aligned_complex


                logger.debug("Successfully align %s onto %s and get the aligned structure %s"%(self._docked_complex, _crystal_receptor, self._aligned_complex))
                #extract the ligand from aligned complex
                try:
                    # extract ligand from output pdb and save to self._aligned_ligand file name
                    # <docked category>-<target id>-<candidate id>_docked_complex_aligned_ligand_<residue name>-<residue id>-<chain id>.pdb
                    extract_ligand_from_complex(self._aligned_complex, self._aligned_ligand, ligand_info = sc_ligand_default)
                    logger.debug("Successfully extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))
                    self._all_aligned_ligand[_crystal_ligand_info] = self._aligned_ligand
                except Exception as ex:
                    logger.exception("Cannot extract ligand file %s from aligned complex %s"%(self._aligned_ligand, self._aligned_complex))

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
                        logger.debug("The rmsd and distance for %s vs %s is %s"%(_crystal_ligand, _docked_ligand_aligned, _rmsd))
                    except Exception as ex:
                        logger.exception("The rmsd calculation for %s vs %s failed"%(_crystal_ligand, _docked_ligand_aligned))
            #get the lowest RMSD and the corresponding distance
            if self._rmsd_dis:
                self._min_rmsd_dis = sorted(self._rmsd_dis.values())[0]


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
            if count % 2 == 0:
                median_index = count / 2
                median = "%-15.3f, " % float(float(sum(sorted_list[median_index-1:median_index+1]))/2)
            else:
                median = "%-15.3f, " % sorted_list[int(round(count / 2))]

        minimum = "%-15.3f, " % min(list_of_value)
        maximum = "%-15.3f, " % max(list_of_value)
    else:
        average = minimum = maximum = median = "%-15s, " % " "

    return average, minimum, maximum, median


class DataContainer(object):

    def __init__(self, ):
        self._data = {}

    def get_all_docked_type(self, docked_type="LMCSS"):
        """Create list of values by examining internal
           data structure of class and for `docked_type` passed
           in return values from all targets registered via
           `register` method.
           :param docked_type: extract values corresponding to this docked_type
           :returns list of values
           """
        docked_type_value = []
        if self._data:
            for target_ID in self._data:
                if docked_type in self._data[target_ID]:
                    docked_type_value.append(self._data[target_ID]
                                             [docked_type])
        return docked_type_value

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

    def layout_json(self, json_filename = "RMSD.json"):
        self._jf = open(json_filename,'w')
        json.dump(self._data, self._jf)
        self._jf.close()

    def layout_plain (self,  plain_filename = "RMSD"):
        self._plain_filename = plain_filename
        self._combined_csv_filename = self._plain_filename + ".csv" 
        self._combined_txt_filename = self._plain_filename + ".txt" 
        self._LMCSS_list = self.get_all_docked_type(docked_type="LMCSS")
        self._SMCSS_list = self.get_all_docked_type(docked_type="SMCSS")
        self._hiResApo_list = self.get_all_docked_type(docked_type="hiResApo")
        self._hiResHolo_list = self.get_all_docked_type(docked_type="hiResHolo")
        self._hiTanimoto_list = self.get_all_docked_type(docked_type="hiTanimoto")

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
    """
    Main entry function
    :param dock_dir:  Path to dock submission.
                      Ex: /..../2018/data...12/stage.6.123_foo.extsubmission
    :param pdb_protein_path: Path to uncompressed pdb
                             Ex: /celpp/pdb.extracted
    :param evaluate_dir: Output directory to use.
                         Ex: /..../stage.7.123_foo.extsubmission.evaluation
    :param blastnfilter_dir: Path to blastnfilter stage
                         Ex: /..../stage.3.blastnfilter
    :param challenge_data_path: Path to challenge data package uncompressed
                                directory
                        Ex: /..../stage.4.challengedata/celpp_week12_2018

    :param update:
    :return:
    """
    #1 copy docked dir info into evaluate dir and create score for targe
    #2 create score dir and copy the ideal files into score dir
    #3 create submit obj and crystal obj in the crystal obj, and store the cyrstal obj directly 
    #4 after each calculation, append the crystal obj to the final result file in the upper level path
    #5 if for individual folder, the crystal obj is empty, also save it and skip it and go to the next case
    #1 switch to evaluation folder
    os.chdir(evaluate_dir)
    current_dir = os.getcwd()

    # DataContainer() is a class to hold the results of the evaluation
    #                  Data is added via register() method and
    #                  the information is held in a dictionary and
    #                  persisted to various files via the layout_json(),
    #                  layout_pickle(), and layout_plain() methods
    result_container = DataContainer()
    result_pickle = "RMSD.pickle"
    pickle_full_path = os.path.join(current_dir, result_pickle) 
    result_json = "RMSD.json"
    json_full_path = os.path.join(current_dir, result_json) 
    result_plain = "RMSD"
    plain_full_path = os.path.join(current_dir, result_plain)

    #get all target dir info from the docked dir
    #
    # This giant for loop iterates through all directories
    # in dock_dir which is stage.6.###_XXX.extsubmission directory
    # and the target_dir is just the directory name. Ex: 1fcz
    target_dirs = list(os.walk(dock_dir))[0][1]
    for target_dir in target_dirs:
        #valid target must have at least one RMSD value return, else this target is invalid
        #
        # This variable is set to True below but never queried or used
        # in the code
        #
        valid_target = False
        
        #copy the target dir to the evaluate dir
        # commands.getoutput("cp -r %s/%s %s"%(dock_dir,target_dir, evaluate_dir))
        #
        # My attempt to replace commands.getoutput with shutil. Looks
        # like the dock_dir/target_dir is copied to the output directory recursively
        # Not sure why since the files could just be referenced.

        target_dir_path = os.path.join(dock_dir, target_dir)
        evaluate_target_dir = os.path.join(evaluate_dir, target_dir)
        logger.debug("Recursively copying " + target_dir_path + " to " + evaluate_target_dir)
        shutil.copytree(target_dir_path, evaluate_target_dir)

        #
        # target_name is totally not needed since target_dir is already
        # the directory name
        #
        target_name = os.path.basename(target_dir)
        logger.info("=================We start to work on this target %s================="%target_name)
        #2 switch into the individual case
        os.chdir(target_name)
        #
        # get_submitted_file_list() looks for all *docked.mol in the directory
        # and generates a tuple (*docked.mol, *docked.pdb)
        #
        all_mol_files, all_receptor_files = get_submitted_file_list()
        if not all_mol_files:
            logger.info("For this folder %s, there is no valid docked structure"%target_dir)
            result_container.register(target_dir, docked_type = None, value = None)
            ######need to register into the result pickle and txt, csv 
            os.chdir(current_dir)
            continue
        else:
            os.mkdir("score")
            # For some reason naming target_name aka target_dir crystal_ID
            #
            #copy the crystal structure first
            crystal_ID = target_name
            #extract the crystal ligand name from blastnfilter result txt file

            # Look at target_dir.txt from stage.3.blastnfilter directory
            # and get the value to right of comma from the line starting with ligand,
            # set that value to "ligand_name" variable below
            #
            blastnfilter_result_name = crystal_ID + ".txt"
            blastnfilter_txt = os.path.join(blastnfilter_dir, blastnfilter_result_name)
            ligand_name = extract_ligand_name(blastnfilter_txt) 
            LMCSS_ligand_name = extract_LMCSS_ligand_name(blastnfilter_txt)

            # copy target_dir aka target_name aka crystal_ID .ent file from pdb
            # directory ie /celpp/pdb.extracted/XX/pdbXXXX.ent
            # and copy to score/crystal.pdb
            #
            #get crystal structure location
            crystal_file = extract_crystal_file(crystal_ID, pdb_protein_path)
            if crystal_file:
                pdbid_local_path = os.getcwd()
                commands.getoutput("cp %s score/crystal.pdb"%(crystal_file))
                try:
                    os.chdir("score")

                    # From score directory call create_crystal_obj which creates a
                    # crystal object. See crystal class above for more information
                    # but in summary creates the following:
                    # target_name/score/crystal.pdb
                    # target_name/score/crystal_ligand-residueId-chainid_ligand.pdb
                    # target_name/score/crystal_ligand-residueId-chainid_receptor.pdb
                    #
                    # Where
                    # crystal.pdb is just the target pdb file obtained from
                    #             pdb database
                    # crystal..._ligand.pdb contains the atoms for the ligand obtained
                    #                       by extracting them out of the crystal.pdb
                    #
                    # crystal..._receport.pdb contains atoms for target and ligand extracted
                    #                         from crystal.pdb
                    #
                    crystal_obj = create_crystal_obj("crystal.pdb", ligand_name, crystal_ID)
                    logger.info("\tSuccessfully create the crystal object")

                    #calculate the distance between the LMCSS ligand and crystal ligand
                    #get the file location

                    # set LMCSS_ori_file_path to <challenge dir>/<target_name>
                    LMCSS_ori_file_path = os.path.join(challenge_data_path, crystal_ID)
                    LMCSS_complex_only_list = []
                    LMCSS_ligand_only_list = []

                    # get a list of all files with LMCSS*.pdb in <challenge dir>/<target name>
                    all_LMCSS_files = glob.glob("%s/LMCSS*.pdb"%LMCSS_ori_file_path)

                    # iterate through all files and if -lig.pdb not in the filename
                    # add it to LMCSS_complex_only_list list otherwise
                    # add it to LMCSS_ligand_only_list
                    #
                    for LMCSS_file in all_LMCSS_files:
                        if "-lig.pdb" not in os.path.basename(LMCSS_file):
                            LMCSS_complex_only_list.append(LMCSS_file)
                        else:
                            LMCSS_ligand_only_list.append(LMCSS_file)

                    # Looks like there should be 1 item in each list otherwise
                    # its an error, but all that happens in an error case is a log message is written
                    # out and code continues
                    #
                    if len(LMCSS_complex_only_list) == 1 and len(LMCSS_ligand_only_list) == 1:
                        LMCSS_complex = LMCSS_complex_only_list[0]
                        LMCSS_ligand = LMCSS_ligand_only_list[0]

                        # below code parses contents of file and returns a string of
                        # <ligand id>-<ligand residue id>
                        # putting it into ligand_info variable
                        #
                        ligand_info = get_ligand_info_from_ligand_file(LMCSS_ligand) 

                        LMCSS_local = os.path.basename(LMCSS_complex)

                        # Copy over the LMCSS*.pdb file that does NOT have -lig.pdb in its name
                        # Usually named: LMCSS-<target>_<target>-<ligand id>.pdb
                        #
                        #copy locally
                        commands.getoutput("cp %s %s"%(LMCSS_complex, LMCSS_local))
                        #create the LMCSS obj
                        # and pass in LMCSS-<target>_<target>-<ligand id>.pdb as first argument
                        # <ligand id>-<ligand residue id> as 2nd argument and the crystal() object
                        # as 3rd argument
                        #
                        LMCSS_obj = LMCSS(LMCSS_local, ligand_info, crystal_obj)
                        #align the
                        LMCSS_obj.align_LMCSS_onto_crystal(check_point_number = 10)
                        #calculate the distance
                        LMCSS_obj.calculate_distance()

                        # The distance between the LMCSS and the crystal is set to LMCSS_distance
                        # variable and registered
                        #
                        LMCSS_distance = LMCSS_obj._min_dis
                        result_container.register(target_dir, docked_type = "LMCSS_ori", value = LMCSS_distance)
                        logger.info("\tSuccessfully calculate the distance between original LMCSS ligand vs crstal ligand. Distance is %s"%LMCSS_distance)
                    else:
                        logger.info("There are %s original LMCSS complex files, which is abnormal and need to check"%len(LMCSS_complex_only_list))
                    os.chdir(pdbid_local_path)
                except Exception as ex:
                    result_container.register(target_dir, docked_type = None, value = None)
                    logger.exception("For this folder %s, could not create the crystal object"%target_dir)
                    os.chdir(current_dir)
                    continue
                #here get the crystal obj and loop all docked structure 
                for mol_index, potential_mol in enumerate (all_mol_files):
                    # all_receptor_files contains all *docked.pdb from the target directory
                    #
                    potential_pdb = all_receptor_files[mol_index]

                    # copy both the .pdb and .mol file to the score directory
                    commands.getoutput("cp %s score"%potential_mol)                     
                    commands.getoutput("cp %s score"%potential_pdb)

                    #change path to local score folder
                    os.chdir("score")
                    local_potential_mol = os.path.basename(potential_mol)               
                    local_potential_pdb = os.path.basename(potential_pdb)

                    # sets docked_structure_type to value of PREFIX in
                    # Files with format as follows:
                    # <PREFIX>-<target id>_<candidate id>_docked.mol
                    # PREFIX should be one of the following according to these docs:
                    # LMCSS, SMCSS, hiResHolo, hiResApo
                    #
                    # See link below for more info
                    # https://github.com/drugdata/D3R/wiki/Challenge-docked-results-file-structure
                    #
                    parsed_name = re.findall('([a-zA-Z0-9]+)-([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_docked.mol',
                                             local_potential_mol)
                    docked_structure_type = parsed_name[0][0]
                    try:
                        # copies local_potential_mol and local_potential_pdb to
                        # score directory and
                        # creates <PREFIX>-<target id>_<candidate id>_docked_ligand.pdb
                        # using schroding structalign
                        docked_obj = docked(local_potential_mol,
                                            local_potential_pdb,
                                            docked_structure_type,
                                            crystal_obj)

                        # creates score/<PREFIX>-<target id>_<candidate id>_docked_ligand_complex.pdb
                        # by merging <PREFIX>-<target id>_<candidate id>_docked_ligand.pdb
                        # with local_potential_pdb using merge_two_pdb() function
                        #
                        docked_obj.create_complex()

                        #
                        # TODO Continue documenting from this point
                        #
                        docked_obj.align_complex_onto_crystal(check_point_number = 10)
                        docked_obj.calculate_rmsd_and_distance()
                        (rmsd, dis) = docked_obj._min_rmsd_dis
                        docked_structure_type_dis = docked_structure_type + "_dis"
                        #rmsd = docked_obj._min_rmsd_dis[0]
                        valid_target = True
                        result_container.register(target_dir, docked_type = docked_structure_type, value = rmsd) 
                        result_container.register(target_dir, docked_type = docked_structure_type_dis, value = dis) 
                        logger.info("\tSuccessfully calculate the rmsd for this category %s, the rmsd is %s"%(docked_structure_type, rmsd))
                        #update the pickle and txt csv file if there is valid case found 
                        result_container.layout_pickle (pickle_filename = pickle_full_path)
                        result_container.layout_json (json_filename = json_full_path)
                        result_container.layout_plain (plain_filename = plain_full_path)
                        os.chdir(pdbid_local_path)
                    except Exception as ex:
                        result_container.register(target_dir, docked_type = docked_structure_type, value = None) 
                        os.chdir(pdbid_local_path)
                        logger.exception("For this type of structure: %s, the rmsd could not be calculated"%docked_structure_type)
                        continue
                logger.info("++++++++++++++++Finished for this target %s+++++++++++++++++"%target_name)
                os.chdir(current_dir)
            else:
                result_container.register(target_dir, docked_type = None, value = None)
                logger.info("For this folder %s, there is no valid crystal structure"%target_dir)
                ######need to register into the result pickle and txt, csv 
                os.chdir(current_dir)
                continue
    #after loop all possible target, layout the pickle and plain files    
    result_container.layout_pickle (pickle_filename = pickle_full_path)
    result_container.layout_json (json_filename = json_full_path)
    result_container.layout_plain (plain_filename = plain_full_path)


def main(args):
    """Main entry into script
    """
    parser = ArgumentParser()
    reqarg = parser.add_argument_group('required named arguments')

    reqarg.add_argument("-d", "--dockdir", metavar="PATH", required=True,
                        help="PATH where we could find the docking "
                             "stage output")

    reqarg.add_argument("-b", "--blastnfilterdir", metavar="PATH",
                        required=True,
                        help="PATH where we could find the blastnfilter "
                             "stage output")

    reqarg.add_argument("-c", "--challengedir", metavar="PATH", required=True,
                        help="PATH where we could find the challengedata "
                             "stage output")

    reqarg.add_argument("-o", "--outdir", metavar="PATH", required=True,
                        help="PATH where we will run the evaluate stage")

    reqarg.add_argument("-p", "--pdbdb", metavar="PATH", required=True,
                        help="PDB DATABANK which we will "
                             "get the crystal structure")

    opt = parser.parse_args(args[1:])

    logging.basicConfig(format='%(asctime)s: %(message)s',
                        datefmt='%m/%d/%y %I:%M:%S',
                        filename=os.path.join(opt.outdir, FINAL_LOG),
                        filemode='w', level=logging.DEBUG)

    logger.debug('Calling main_score function')
    main_score(opt.dockdir, opt.pdbdb, opt.outdir,
               opt.blastnfilterdir, opt.challengedir)
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))

