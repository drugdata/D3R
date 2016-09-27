#!/usr/bin/env python

import commands
import os
import glob
import logging
import time
from openeye.oechem import *
import pickle
import numpy
import re

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )


schrodingerScriptText = '''
import schrodinger.structutils.transform as ssut
import schrodinger.structure as ss
import numpy as np
import sys
if len(sys.argv) != 4:
    print "Usage: python temp.py inputmol outputmol rotationmatrix"
inFile = sys.argv[1]
outFile = sys.argv[2]
matrixFile = sys.argv[3]
rotMat = np.genfromtxt(matrixFile, skip_header=1)
#rotMat = np.concatenate((rotMat.T, np.zeros((1,4))))
rotMat = np.concatenate((rotMat, np.zeros((4,1))), axis=1)
rotMat[0,3] = rotMat[3,0]
rotMat[1,3] = rotMat[3,1]
rotMat[2,3] = rotMat[3,2]
#print rotMat
mol = ss.StructureReader(inFile).next()
ssut.transform_structure(mol, rotMat)
writer = ss.StructureWriter(outFile)
writer.append(mol)
writer.close()
'''


def merge_two_pdb (receptor, ligand, complex_pdb):
    complex_lines = []
    f1 = open(receptor, "r")
    protein_lines = f1.readlines()
    f1.close()
    for p_line in protein_lines:
        if p_line [:6] not in ["CONECT", "ENDMDL", "END   "]:
            complex_lines.append(p_line)
    complex_lines.append("TER   \n")
    f2 = open(ligand, "r")
    ligand_lines = f2.readlines()
    f2.close()
    for l_line in ligand_lines:
        if l_line [:6] not in ["REMARK", "MODEL ", "CONECT", "ENDMDL", "END   "]:
            complex_lines.append(l_line) 
    f3 = open(complex_pdb, "w")
    f3.writelines(complex_lines)
    f3.close()

def heavy_atom (mol):
#check the heavy atom numbers
    heavy_atom_numbers = 0
    for atom in mol.GetAtoms():
        if not atom.IsHydrogen():
            heavy_atom_numbers += 1
    return heavy_atom_numbers

def mcs(ref_mol, fit_mol, ref_ignore = [], fit_ignore = []):
#do the mcs search and return OEMatch object.
    #count heavy atoms 
    ref_mol_heavy_atoms = heavy_atom(ref_mol)
    fit_mol_heavy_atoms = heavy_atom(fit_mol)
    #ignore Hydrogen 
    OESuppressHydrogens(fit_mol)
    OESuppressHydrogens(ref_mol)
    #set atom and bond expression                                       
    atomexpr = OEExprOpts_AtomicNumber                                  
    bondexpr = 0                                                        
    #do the mcs search, using defined atom, bond expression options, and make sure using the small molecule as reference mol
    if not ref_mol_heavy_atoms > fit_mol_heavy_atoms:                   
        mcss = OEMCSSearch( ref_mol, atomexpr, bondexpr, True)          
    else:                                                               
        mcss = OEMCSSearch( fit_mol, atomexpr, bondexpr, True)          
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
        match1_heavy_atom = heavy_atom(mol1)                            
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
            match2_heavy_atom = heavy_atom(mol2)
            ofs2 = oemolostream("match2_%s.pdb"%j)                      
            OEWriteMolecule(ofs2, mol2)                                 
            ofs2.close()                                                
            if match1_heavy_atom != match2_heavy_atom:                  
                print "WARNING: match1_%s.pdb and match2_%s.pdb have different heavy atoms, need to check..."%(i,j)
            for mp1, mp2 in zip(match1.GetAtoms(), match2.GetAtoms()):  
                ref_name = mp1.target.GetName().strip()                 
                fit_name = mp2.target.GetName().strip()                 
                if (ref_name not in ref_ignore) and (fit_name not in fit_ignore ) :   
                    new_match.AddPair (mp1.target, mp2.target)          
                else:                                                   
                    print "             Deleted atom pair", mp1.target.GetName(), mp2.target.GetName()
            new_match_list.append(new_match)                            
            new_match_dic[new_match] = (["match1_%s_vs_match2_%s"%(i,j), check_list ])
    return new_match_dic 

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def add_symbol (pdb_file):
    # check the file format, if it is pdb file then do the add.
    if os.path.splitext(pdb_file)[-1] == ".pdb":
        f = open(pdb_file, "r")
        lines = f.readlines()
        new_lines = []
        for line in lines:
            #lines with atom informations  
            if "ATOM" or "HETATM" in line.split()[0]:
                if is_number(line.split()[-1]) and is_number(line.split()[-2]):
                    #ensure the end has no atom symbol, which means the end is number. Also ensure the end is not charge. 
                    atom_name = line.split()[2]
                    #remove numbers in atom name.                       
                    atom_symbol = ''.join(i for i in atom_name if not i.isdigit())     
                    #add whitespace because of the requirment of the pdb format       
                    lenth_symbol = "          " + "%s"%atom_symbol      
                    new_line = line.split("\n")[0].strip()+ lenth_symbol + "\n"       
                    new_lines.append(new_line)                          
                else:                                                   
                    new_lines.append(line)                              
            else:                                                       
                new_lines.append(line)                                  
        file_name, file_extension = os.path.splitext(pdb_file)          
        out_file = file_name + "_new_atom_symbol" + file_extension      
        f2 = open(out_file,"w")                                         
        f2.writelines(new_lines)
        return out_file                                                 
    else:                                                               
        #if the format is not pdb then return original                  
        return pdb_file

def main_rmsd(ref_struc, fit_struc):
    reffs = oemolistream()
    reffs.open(add_symbol(ref_struc))
    fitfs = oemolistream() 
    fitfs.open(add_symbol(fit_struc))
    refmol = OEGraphMol()
    OEReadMolecule(reffs, refmol)
    for fitmol in fitfs.GetOEGraphMols():
        rmsd  = OERMSD(refmol, fitmol)
        rot  = OEDoubleArray(9)
        trans = OEDoubleArray(3)
        automorph = True
        heavyOnly = True
        overlay = True
        rmsd  = OERMSD(refmol, fitmol, automorph, heavyOnly, overlay, rot, trans)
        ss = mcs(refmol, fitmol,)
        mcss_rmsd_list = []
        match_info = []
        for mcss in ss.keys():
            #mcss_rmsd = OERMSD(refmol, fitmol, mcss, overlay, rot, trans)
            mcss_rmsd = OERMSD(refmol, fitmol, mcss)
            mcss_rmsd_list.append(mcss_rmsd)
            match_info.append(ss[mcss])
        for filename in glob.glob('*_new_atom_symbol.pdb'):             
            os.remove(filename) 
        SMCSS_number_index = mcss_rmsd_list.index(min(mcss_rmsd_list))
        biggest_number_index = mcss_rmsd_list.index(max(mcss_rmsd_list))
    return min(mcss_rmsd_list)


def old_structure_align(template_pro, target_pro):
    commands.getoutput("$SCHRODINGER/utilities/structalign %s %s"%(template_pro, target_pro))
    out_rot_target_pro = "rot-" + target_pro
    if os.path.isfile(out_rot_target_pro):
        return out_rot_target_pro
    else:
        logging.info("Could not algin %s onto %s"%(target_pro,template_pro))

def structure_align(prefix, actual_xtal_pdb, receptor_in, ligand_in):
    receptor_mae = 'rot_%s_receptor1.mae' %(prefix)
    receptor_out = 'rot_%s' %(receptor_in)
    ligand_out = 'rot_%s' %(ligand_in)
    #matrixFile = receptor_in.replace('.pdb','.rot')
    matrixFile = 'temp_receptor.rot'
    with open('temp.py','wb') as of:
        of.write(schrodingerScriptText)
    
    
    ## On the old version of maestro/structalign that's on nif1, we can't get the rotation matrix by aligning two pdbs. Until it gets updated, I'll convert a copy to mae and get the rotation matrix by aligning that.
    
    cmd0 = '$SCHRODINGER/utilities/structconvert -ipdb %s -omae temp_xtal.mae >& %s-0-convertXtalPdbToMae' %(actual_xtal_pdb, prefix)
    commands.getoutput(cmd0)
    cmd0point5 = '$SCHRODINGER/utilities/structconvert -ipdb %s -omae temp_receptor.mae >& %s-0.5-convertCanPdbToMae' %(receptor_in, prefix)
    commands.getoutput(cmd0point5)
    
    ## NOTE! the default rotation matrix output by structalign only has 3 decimal places. This might make a difference. 
    ## To fix it, make a copy of $SCHRODINGER/mmshare-v<whatever>/bin/Linux-x86_64/structalign_utility.py in which the %.3f's in the mat.write lines are replaced with %.5f's.
    #cmd1 = '$SCHRODINGER/utilities/structalign -matrix %s %s >& %s_structAlignOut' %(actual_xtal_pdb, receptor_in, prefix)
    cmd1 = '$SCHRODINGER/utilities/structalign -matrix temp_xtal.mae temp_receptor.mae >& %s-1-alignMatrixGenOut' %(prefix)
    commands.getoutput(cmd1)
    
    cmd2 = '$SCHRODINGER/utilities/python temp.py %s %s %s >& %s-2-receptorAlignOut' %(receptor_in, receptor_mae, matrixFile, prefix)
    commands.getoutput(cmd2)
    
    cmd3 = '$SCHRODINGER/utilities/python temp.py %s %s %s >& %s-3-ligandAlignOut' %(ligand_in, ligand_out, matrixFile, prefix)
    commands.getoutput(cmd3)
    
    cmd4 = '$SCHRODINGER/utilities/structconvert -imae %s -opdb %s >& %s-4-structConvertOut' %(receptor_mae, receptor_out, prefix)
    commands.getoutput(cmd4)
    
    return receptor_out, ligand_out

def make_complex_pdb(receptor_pdb, ligand_mol, complex_pdb):
    ## First, convert the ligand to pdb
    ligand_pdb = ligand_mol.replace('.mol','_lig.pdb')
    commands.getoutput('babel -imol %s -opdb %s' %(ligand_mol, ligand_pdb))
    ## Now combine the ligand and receptor pdbs
    commands.getoutput('babel --join -ipdb %s -ipdb %s -opdb %s' %(receptor_pdb, ligand_pdb, complex_pdb))
    return ligand_pdb, complex_pdb
 

def layout_result (pickle_file, out_file, outformat = "txt"):
    data = []
    p_file = open(pickle_file,"r")
    score_dic = pickle.load(p_file)
    p_file.close()
    if outformat == "csv":
        data = ["%-20s%-10s, %-10s, %-10s, %-10s, %-10s \n"%("Target_PDBID,", "LMCSS", "SMCSS", "hiResApo", "hiResHolo", "hiTanimoto")]
    elif outformat =="txt":
        data = ["%-20s%-10s%-10s%-10s%-10s%-10s \n"%("Target_PDBID", "LMCSS", "SMCSS", "hiResApo", "hiResHolo", "hiTanimoto")]
    #add hiTanimoto 0909 sliu
    all_pro_type = ["LMCSS", "SMCSS", "hiResApo", "hiResHolo", "hiTanimoto"]
    LMCSS_list =  []
    SMCSS_list = []
    hiResApo_list = []
    hiResHolo_list = []
    hiTanimoto_list = []
    abnormal_dic = {}
    total_pdb = 0
    for pdbid in score_dic:
        for pro_type in all_pro_type:
            if pro_type in score_dic[pdbid]:
                if pro_type == "LMCSS":
                    LMCSS_list.append(score_dic[pdbid][pro_type])
                    if float(score_dic[pdbid][pro_type]) > 8.0:
                        abnormal_dic[pdbid] = score_dic[pdbid][pro_type]
                if pro_type == "SMCSS":
                    SMCSS_list.append(score_dic[pdbid][pro_type])
                if pro_type == "hiResApo":
                    hiResApo_list.append(score_dic[pdbid][pro_type]) 
                if pro_type == "hiResHolo":
                    hiResHolo_list.append(score_dic[pdbid][pro_type])
                if pro_type == "hiTanimoto":
                    hiTanimoto_list.append(score_dic[pdbid][pro_type])
        total_pdb += 1
    if outformat == "csv":
        data.append("%-20s%-10s, %-10s, %-10s, %-10s, %-10s, \n"%("Number_of_cases,", len(LMCSS_list), len(SMCSS_list), len(hiResApo_list), len(hiResHolo_list), len(hiTanimoto_list)))
    elif outformat == "txt":
        data.append("%-20s%-10s%-10s%-10s%-10s%-10s\n"%("Number_of_cases", len(LMCSS_list), len(SMCSS_list), len(hiResApo_list), len(hiResHolo_list), len(hiTanimoto_list)))
    
    whole_list = (LMCSS_list, SMCSS_list, hiResApo_list, hiResHolo_list, hiTanimoto_list)
    if outformat == "csv":
        average_line, max_line, min_line = ("%-20s"%"Average,", "%-20s"%"Maximum,", "%-20s"%"Minimum,")
    elif outformat == "txt":
        average_line, max_line, min_line = ("%-20s"%"Average", "%-20s"%"Maximum", "%-20s"%"Minimum")
    for index, this_list in enumerate (whole_list):
        if len(this_list) == 0:
            if outformat == "csv":
                average_line += "%-10s, "%(" ")
                min_line += "%-10s, "%(" ")
                max_line += "%-10s, "%(" ")
            elif outformat == "txt":
                average_line += "%-10s"%(" ")
                min_line += "%-10s"%(" ")
                max_line += "%-10s"%(" ")
        else:
            if outformat == "csv":
                average_line += "%-10.3f, "%(numpy.average(this_list))
                min_line += "%-10.3f, "%(min(this_list))
                max_line += "%-10.3f, "%(max(this_list))
            if outformat == "txt":
                average_line += "%-10.3f"%(numpy.average(this_list))
                min_line += "%-10.3f"%(min(this_list))
                max_line += "%-10.3f"%(max(this_list))
    average_line += "\n" 
    max_line += "\n"
    min_line += "\n"
    data.append(average_line)
    data.append(min_line)
    data.append(max_line)
       
    #add main score lines 
    for pdbid in score_dic:
        new_line_score = ""
        valid_line = False
        for pro_type in all_pro_type:
            if pro_type in score_dic[pdbid]:
                if outformat == "csv":
                    new_line_score += "%-10.3f, "%score_dic[pdbid][pro_type] 
                if outformat == "txt":
                    new_line_score += "%-10.3f"%score_dic[pdbid][pro_type] 
                valid_line = True
            else:
                if outformat == "csv":
                    new_line_score += "%-10s, "%(" ")
                if outformat == "txt":
                    new_line_score += "%-10s"%(" ")
        if valid_line:
            if outformat == "csv":
                pdbid_full = "%s,"%pdbid
            if outformat == "txt":
                pdbid_full = "%s"%pdbid
            new_line = "%-20s"%(pdbid_full)
            new_line += new_line_score
            new_line += "\n" 
            data.append(new_line)
    #data.append("=====Total number of target protein: %s =====\n"%total_pdb)
    #append total number of valid pdb for each type 
    #data.append("%-20s%-10s, %-10s, %-10s, %-10s, %-10s, \n"%("Valid_cases,", len(LMCSS_list), len(SMCSS_list), len(hiResApo_list), len(hiResHolo_list), len(hiTanimoto_list)))
    #data.append("%-20s%-10.3f, %-10.3f, %-10.3f, %-10.3f, %-10.3f, \n"%("Average,", numpy.average(LMCSS_list), numpy.average(SMCSS_list), numpy.average(hiResApo_list), numpy.average(hiResHolo_list), numpy.average(hiTanimoto_list)))
    #data.append("=====Abnormal RMSDs (LMCSS cases where the RMSD > 8.0)=====\n")
    #for abnormal_id in abnormal_dic:
    #    abnormal_id_full = "\"%s\","%abnormal_id
    #    data.append("%-20s%-10.3f, \n"%(abnormal_id_full, abnormal_dic[abnormal_id]))
    out_txt = open(out_file, "w")
    out_txt.writelines(data)
    out_txt.close()  
     

def main_score (dock_dir, pdb_protein_path, evaluate_dir, update= True):
    pickle_file_name = "RMSD.pickle"
    pickle_out = open (pickle_file_name, "w")
    score_dic = {}
    os.chdir(evaluate_dir)
    current_dir = os.getcwd()
    #target_dirs = []
    #all_pdbids = next(os.walk(dock_dir))[1]    
    ##all_pdbids = ['4xm7']
    #for all_pdb_path in os.walk(dock_dir):
    #    if all_pdb_path[0] != dock_dir:
    #        
    #        if all_pdb_path[0].split("/")[-1] in all_pdbids:
    #            #make sure that we didn't visit the deep folders other than the pdbid ones
    #            target_dirs.append(all_pdb_path[0])
    
    ## os.walk will return a tuple of (directory absolute path, subdirectories, files)

    target_dirs = list(os.walk(dock_dir))[0][1]

    
    for target_dir in target_dirs:
        all_docked_structures = []
        valid_target = False # this flag indicates whether we have at least one valid docked structure for scoring
        commands.getoutput("cp -r %s/%s %s"%(dock_dir,target_dir, evaluate_dir))
        target_name = os.path.basename(target_dir)
        os.chdir(target_name)
        ## Ensure the target is in score_dic
        if target_name not in score_dic:
            score_dic[target_name] = {}
            
        #now we are in folders named by pdbid
        current_dir_layer_2 = os.getcwd()
        ##################
        #Do the scoring here
        #1, get all the docked structures and crystal structure
        potential_mols = glob.glob('*/*docked.mol' )
        ## Go copy in all the submitted poses
        for potential_mol in potential_mols:
            
            ## First, ensure that this ligand has a corresponding receptor
            potential_receptor = potential_mol.replace('.mol','.pdb')
            if not(os.path.isfile(potential_receptor)):
                logging.info('There is no receptor corresponding to docked molecule %s. Skipping scoring this pose.' %(potential_mol))
                continue
            
            ## If it does, set up the directory for scoring 
            if valid_target == False:
                valid_target = True
                os.mkdir("score")
                
            ## don't copy in the original pdb and mol files - Just convert them to mae
            #complex_mae = os.path.basename(potential_mol).replace('.mol','_complex.mae')
            #commands.getoutput('$SCHRODINGER/utilities/structcat -imol score/%s -ipdb score/%s -omae score/%s' %(potential_mol, potential_receptor, complex_mae))
            #all_docked_structures.append(complex_mae)
            
            ## copy in the structures intended for scoring
            commands.getoutput("cp %s score" %potential_mol)
            commands.getoutput("cp %s score" %potential_receptor)
            
            ## Combine them into a single structure for RSR calculation
            ## The score subdirectory is flat, so we need to drop the target directory structure attached to potential_mol
            local_potential_mol = os.path.basename(potential_mol)
            local_potential_receptor = os.path.basename(potential_receptor)
            #intermediate_lig_pdb = local_potential_mol.replace('.mol','_ligand.pdb')
            ## The complex mae will be called "<method>-<target>_<candidate>_docked_complex.pdb"
            #complex_pdb = local_potential_mol.replace('.mol','_complex.pdb')
            #commands.getoutput('babel -imol score/%s -opdb score/%s' %(potential_mol, intermediate_lig_pdb))
            #commands.getoutput('babel --join -ipdb score/%s -ipdb score/%s -opdb score/%s' %(potential_receptor, intermediate_lig_pdb, complex_pdb))
            
            
            ## Package the ligand, receptor, and complex file up for the next step
            all_docked_structures.append((local_potential_mol, local_potential_receptor))
        
        ## Check whether it makes sense to continue scoring
        if not(valid_target):
            #no docked structure available for this case, skip this one
            logging.info("There are no valid docked structures for this PDBID: %s"%target_name)
            os.chdir(current_dir)
            continue
        
        ## Extract the crystal (local pdb database is organized by middle 2 characters of pdb codes)
        crystal_id = target_name
        crystal_pdb_folder_name = crystal_id[1:3]    
        crystal_ent_file = "pdb" + crystal_id + ".ent"
        crystal_pdbloc = os.path.join(pdb_protein_path,crystal_pdb_folder_name, crystal_ent_file )
        ## Ensure that the xtal entry for this pdb code exists
        if not os.path.isfile(crystal_pdbloc):
            logging.info("Unable to find the ent file associate with the crystal pdb: %s at location %s"%(crystal_id, crystal_pdbloc))
            os.chdir(current_dir)
            continue
        else:
            commands.getoutput("cp %s score/crystal.pdb"%(crystal_pdbloc))
            
        ## Finish step 1 
        ## Move into score dir 
        os.chdir("score")    
        logging.info("=============Start working in this case:%s=========="%target_name)
        #2, split the maegz file into pdb files and combine receptor with ligand
        #align the complex onto the crystal structure and split into ligands with receptor
        
        #all_docked_structures = glob.glob("*_docked.mol") 
        ## remove all files start with rot-
        #previous_files = glob.glob("rot-*_docked.mol")
        #for item in previous_files:
        #    all_docked_structures.remove (item)
        
        ## Split the xtal ligands from the xtal stucture
        commands.getoutput("$SCHRODINGER/run split_structure.py -many_files -m ligand crystal.pdb crystal.pdb")
        if not os.path.isfile("crystal_ligand1.pdb"):
            logging.info("Crystal structure is not splittable...")
            os.chdir(current_dir_layer_2)
            os.chdir(current_dir)
            continue    
        else:
            # Get the filenames of all the crystal ligands
            crystal_ligand_list = glob.glob("crystal_ligand*.pdb")
        
        ## Remember that all_docked_structures is a list of tuples: (potential_mol, potential_receptor, complex_pdb))
        for docked_lig_mol, docked_receptor in all_docked_structures:
            
            ## Get the method type, target, and candidate info from the filename
            # for example, this will parse 'hiResApo-5hib_2eb2_docked.mol' into [('hiResApo', '5hib', '2eb2')]
            parsed_name = re.findall('([a-zA-Z0-9]+)-([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_docked.mol',docked_lig_mol)
            if len(parsed_name) != 1:
                logging.info('Failed to parse docked structure name "%s". Parsing yielded %r' %(docked_lig_mol, parsed_name))
                continue
            docked_structure_type = parsed_name[0][0]
            docked_structure_target = parsed_name[0][1]
            docked_structure_candidate = parsed_name[0][2]
            
            
            ## Align the submission to the crystal structure
            file_prefix = "%s-%s_%s_docked" %(docked_structure_type,
                                              docked_structure_target,
                                              docked_structure_candidate)             
            aln_complex_pdb = '%s_complex.pdb' %(file_prefix)
            
            aln_recep, aln_lig = structure_align(file_prefix, 'crystal_receptor1.pdb', docked_receptor, docked_lig_mol)
            aln_lig_pdb, aln_complex_pdb = make_complex_pdb(aln_recep, aln_lig, aln_complex_pdb)
            #docked_lig_pdb = docked_lig_mol.replace('.mol','.pdb')
            
            
            ## Now compare the RMSDs of every ligand to those in the crystal
            ## We need the ligand as a pdb
            try:
                rmsd_list = []
                for crystal_ligand in crystal_ligand_list:
                    rmsd = main_rmsd(crystal_ligand, aln_lig_pdb)
                    logging.info( "RMSD for the first ligand: %s comparing with crystal ligand :%s, is : %s"%(aln_lig_pdb, crystal_ligand, rmsd))
                    rmsd_list.append(rmsd)
                if docked_structure_type not in score_dic[target_name]:
                    score_dic[target_name][docked_structure_type] = min(rmsd_list)
            except:
                logging.info("RMSD cannot be calculated for the ligand: %s"%(aln_lig_pdb))
            #structure_type = all_docked_structure.split("_")[0]
            
            #first split, then combine together, then do the aligned
            #all_docked_structure_title = all_docked_structure.split(".")[0]
            #all_docked_structure_pdb = all_docked_structure_title + ".pdb"
            #commands.getoutput("$SCHRODINGER/run split_structure.py -many_files -m ligand %s %s"%(all_docked_structure, all_docked_structure_pdb))
            
            # should get LMCSS_dock_pv_receptor1.pdb and LMCSS_dock_pv_ligand1.pdb, may get LMCSS_dock_pv_ligand2.pdb
            # now combine receptor and ligand
            #all_docked_structure_receptor = all_docked_structure_title + "_receptor1.pdb"
            
            #here, just focus on the first ligand now, probably will come back to get all other ligands sliu 2/12/2016
            #all_docked_structure_ligand = all_docked_structure_title + "_ligand1.pdb"
            #all_docked_structure_complex = all_docked_structure_title + "_complex1.pdb"
            
            #here once split, there may not have the ligand structure, because the ligand is big and schrodinger regard this ligand as a receptor, so no ligand detected... need to come back this this issue later sliu 2/25/2016, tmp fix, add an exception 
            #try:
            #    merge_two_pdb(all_docked_structure_receptor, all_docked_structure_ligand, all_docked_structure_complex)
            #except:
            #    logging.info("Could not merge %s and %s together, need to check"%(all_docked_structure_ligand, all_docked_structure_complex))
            #if not os.path.isfile(all_docked_structure_complex):
            #    logging.info("Could not get combined complex file:%s"%all_docked_structure_complex)
            #    continue
            #aligned_sturture = structure_align("crystal.pdb", all_docked_structure_complex)
            #if aligned_sturture:
            #    commands.getoutput("$SCHRODINGER/run split_structure.py -many_files -m ligand %s %s"%(aligned_sturture,aligned_sturture))
            #    aligned_title = aligned_sturture.split(".")[0]
            #    #get the splitted ligand and protein
            #    ligand_list = glob.glob("%s_ligand*.pdb"%aligned_title)
            #    #just use the first ligand here
            #    top_ligand_pdb = "%s_ligand1.pdb"%aligned_title
            #    try:
            #        #change to the rmsd list style, becaue the crystal structure may have multiple ligand
            #            #choose the lowest rmsd to store.
            #        rmsd_list = []
            #        for crystal_ligand in crystal_ligand_list:
            #            rmsd = main_rmsd("%s"%crystal_ligand, "%s"%top_ligand_pdb)        
            #            rmsd_list.append(rmsd)
            #            logging.info( "RMSD for the first ligand: %s comparing with crystal ligand :%s, is : %s"%(top_ligand_pdb, crystal_ligand, rmsd))
            #        if structure_type not in score_dic[scoreable_path_local]:
            #            score_dic[scoreable_path_local][structure_type] = min(rmsd_list)
            #    except:
            #        logging.info("RMSD cannot be calculate for the ligand: %s"%top_ligand_pdb)
        #Finish scoring, come back to the folder with pdbid
        os.chdir(current_dir_layer_2)
        # Finish scoring for this pdbid, come back to the main folder
        os.chdir(current_dir)
        ##################
    pickle.dump(score_dic, pickle_out)
    pickle_out.close()
    return pickle_file_name

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-d", "--dockdir", metavar="PATH",
                      help="PATH where we could find the docking stage output")

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
    # running under this dir
    running_dir = os.getcwd()
    pickle_result = main_score(dockDir, pdbloc, evaluateDir, )
    pickle_loc = os.path.join(running_dir, "RMSD.pickle")
    txt_result = layout_result(pickle_loc, "RMSD.txt", outformat = "txt")
    txt_result = layout_result(pickle_loc, "RMSD.csv", outformat = "csv")
    # move the final log file to the result dir
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s" % (log_file_path, evaluateDir))
    pickle_file_path = os.path.join(running_dir, 'RMSD.pickle')
    commands.getoutput("mv %s %s" % (pickle_file_path, evaluateDir))
