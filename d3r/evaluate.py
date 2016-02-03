#!/usr/bin/evn python

import commands
import os
import glob
import logging
import time
from openeye.oechem import *
import pickle
import numpy

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )

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
                print "WARNING: match1_%s.pdb and match2_%s.pdb have different heavy atoms, need to check..."
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
        smallest_number_index = mcss_rmsd_list.index(min(mcss_rmsd_list))
        biggest_number_index = mcss_rmsd_list.index(max(mcss_rmsd_list))
    return min(mcss_rmsd_list)


def structure_align(template_pro, target_pro):
    commands.getoutput("$SCHRODINGER/utilities/structalign %s %s"%(template_pro, target_pro))
    out_rot_target_pro = "rot-" + target_pro
    if os.path.isfile(out_rot_target_pro):
        return out_rot_target_pro
    else:
        logging.info("Could not algin %s onto %s"%(target_pro,template_pro))


def layout_result (pickle_file, txt_file):
    data = []
    p_file = open(pickle_file,"r")
    score_dic = pickle.load(p_file)
    p_file.close()
    data = ["%-20s%-10s%-10s%-10s%-10s \n"%(" ", "largest", "smallest", "apo", "holo")]
    all_pro_type = ["largest", "smallest", "apo", "holo"]
    largest_list =  []
    smallest_list = []
    apo_list = []
    holo_list = []
    abnormal_dic = {}
    total_pdb = 0
    for pdbid in score_dic:
        for pro_type in all_pro_type:
            if pro_type in score_dic[pdbid]:
                if pro_type == "largest":
                    largest_list.append(score_dic[pdbid][pro_type])
                    if float(score_dic[pdbid][pro_type]) > 8.0:
                        abnormal_dic[pdbid] = score_dic[pdbid][pro_type]
                if pro_type == "smallest":
                    smallest_list.append(score_dic[pdbid][pro_type])
                if pro_type == "apo":
                    apo_list.append(score_dic[pdbid][pro_type]) 
                if pro_type == "holo":
                    holo_list.append(score_dic[pdbid][pro_type])
        total_pdb += 1
    for pdbid in score_dic:
        new_line_score = ""
        for pro_type in all_pro_type:
            if pro_type in score_dic[pdbid]:
                new_line_score += "%-10.3f"%score_dic[pdbid][pro_type] 
        if new_line_score:
            new_line = "%-20s"%pdbid
            new_line += new_line_score
            new_line += "\n" 
            data.append(new_line)
    data.append("=====Total number of query protein: %s =====\n"%total_pdb)
    #append total number of valid pdb for each type 
    data.append("%-20s%-10s%-10s%-10s%-10s\n"%("Valid cases", len(largest_list), len(smallest_list), len(apo_list), len(holo_list)))
    data.append("%-20s%-10.3f%-10.3f%-10.3f%-10.3f\n"%("Average", numpy.average(largest_list), numpy.average(smallest_list), numpy.average(apo_list), numpy.average(holo_list)))
    data.append("=====Abnormal RMSDs =====\n")
    for abnormal_id in abnormal_dic:
        data.append("%-20s%-10.3f\n"%(abnormal_id, abnormal_dic[abnormal_id]))
    out_txt = open(txt_file, "w")
    out_txt.writelines(data)
    out_txt.close()  
     

def main_score (stage_4_result, pdb_protein_path, stage_5_working, update= True):
    pickle_file_name = "RMSD.pickle"
    pickle_out = open (pickle_file_name, "w")
    score_dic = {}
    os.chdir(stage_5_working)
    current_dir = os.getcwd()
    scoreable_paths = []
    all_pdbids = next(os.walk(stage_4_result))[1]    
    #all_pdbids = ['4xm7']
    for all_pdb_path in os.walk(stage_4_result):
        if all_pdb_path[0] != stage_4_result:
            
            if all_pdb_path[0].split("/")[-1] in all_pdbids:
                #make sure that we didn't visit the deep folders other than the pdbid ones
                scoreable_paths.append(all_pdb_path[0])
    for scoreable_path in scoreable_paths:
        commands.getoutput("cp -r %s %s"%(scoreable_path, stage_5_working))
        scoreable_path_local = os.path.basename(scoreable_path)
        os.chdir(scoreable_path_local)
        if scoreable_path_local not in score_dic:
            score_dic[scoreable_path_local] = {}
        print "SSSSSSSSSSSSS", score_dic
        #now we are in folders named by pdbid
        current_dir_layer_2 = os.getcwd()
        ##################
        #Do the scoring here
        #1, get all the docked strutures and crystal structure 
        potential_structures = []
        for titles in ["largest", "smallest", "apo", "holo"]:
            potential_structures.append("%s/%s_dock_pv.maegz"%(titles,titles))
        for potential_structure in potential_structures:
            if os.path.isfile(potential_structure):
                #Found one docked struture, trigue the scoring
                if not os.path.isdir("score"):
                    os.mkdir("score")
                commands.getoutput("cp %s score"%potential_structure)
        if not os.path.isdir("score"):
            #no docked structue available for this case, skip this one
            logging.info("There is no docked struture for this PDBID: %s"%scoreable_path_local)
            os.chdir(current_dir)
            continue
        #extract the crystal 
        crystal_id = scoreable_path_local
        crystal_pdb_folder_name = crystal_id[1:3]    
        crystal_ent_file = "pdb" + crystal_id + ".ent"
        crystal_pdbloc = os.path.join(pdb_protein_path,crystal_pdb_folder_name, crystal_ent_file )
        if not os.path.isfile(crystal_pdbloc):
            logging.info("Unable to find the ent file associate with the crystal pdb: %s at location %s"%(crystal_id, crystal_pdbloc))
            os.chdir(current_dir)
            continue
        else:
            commands.getoutput("cp %s score/crystal.pdb"%(crystal_pdbloc))
        #finish step 1 
        #go to score dic and do score 
        os.chdir("score")    
        logging.info("=============Start working in this case:%s=========="%scoreable_path_local)
        #2, align all maegz file onto the crystal structure and split right after alignment
        #get all the protein and align the protein to the crystal structure
        all_docked_structures = glob.glob("*.maegz") 
        previous_files = glob.glob("rot-*.maegz")
        for item in previous_files:
            all_docked_structures.remove (item)
        #remove all files start with rot-
        #get the template ligand
        commands.getoutput("$SCHRODINGER/run split_structure.py -many_files -m ligand crystal.pdb crystal.pdb")
        if not os.path.isfile("crystal_ligand1.pdb"):
            logging.info("Crystal structure is not splitable...")
            os.chdir(current_dir_layer_2)
            os.chdir(current_dir)
            continue    
        for all_docked_structure in all_docked_structures:
            structure_type = all_docked_structure.split("_")[0]
            aligned_sturture = structure_align("crystal.pdb", all_docked_structure)
            if aligned_sturture:
                commands.getoutput("$SCHRODINGER/run split_structure.py -many_files -m ligand %s %s"%(aligned_sturture,aligned_sturture))
                aligned_title = aligned_sturture.split(".")[0]
                #get the splitted ligand and protein
                ligand_list = glob.glob("%s_ligand*.maegz"%aligned_title)
                #just use the first ligand here
                top_ligand = "%s_ligand1.maegz"%aligned_title
                top_ligand_pdb = top_ligand.split(".")[0]+".pdb"
                #convert to pdb file
                commands.getoutput("$SCHRODINGER/utilities/pdbconvert -imae %s -opdb %s"%(top_ligand, top_ligand_pdb)) 
                try:
                    rmsd = main_rmsd("crystal_ligand1.pdb", "%s"%top_ligand_pdb)        
                    logging.info( "RMSD for the first ligand: %s is : %s"%(top_ligand, rmsd))
                    if structure_type not in score_dic[scoreable_path_local]:
                        score_dic[scoreable_path_local][structure_type] = rmsd
                except:
                    logging.info("RMSD cannot be calculate for the ligand: %s"%top_ligand)
        # Finsih scoring, come back to the folder with pdbid
        os.chdir(current_dir_layer_2)
        # Finsih scoring for this pdbid, come back to the main folder
        os.chdir(current_dir)
        ##################
    pickle.dump(score_dic, pickle_out)
    pickle_out.close()
    return pickle_file_name

if ("__main__") == (__name__):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-d", "--dockdir", metavar="PATH",
                      help="PATH where we could find the stage 4 docking "
                           "output")

    parser.add_option("-o", "--outdir", metavar="PATH",
                      help="PATH where we run stage 5")

    parser.add_option("-p", "--pdbdb", metavar="PATH",
                      help="PDB DATABANK which we will "
                           "get the crystal structure")
    logger = logging.getLogger()
    logging.basicConfig(format='%(asctime)s: %(message)s',
                        datefmt='%m/%d/%y %I:%M:%S', filename='final.log',
                        filemode='w', level=logging.INFO)
    (opt, args) = parser.parse_args()
    stage_4_result = opt.dockdir
    stage_5_result = opt.outdir
    pdbloc = opt.pdbdb
    # running under this dir
    running_dir = os.getcwd()
    pickle_result = main_score(stage_4_result, pdbloc, stage_5_result, )
    pickle_loc = os.path.join(running_dir, "RMSD.pickle")
    txt_result = layout_result(pickle_loc, "RMSD.txt")
    # move the final log file to the result dir
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s" % (log_file_path, stage_5_result))
    pickle_file_path = os.path.join(running_dir, 'RMSD.pickle')
    commands.getoutput("mv %s %s" % (pickle_file_path, stage_5_result))
