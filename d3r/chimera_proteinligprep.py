#!/usr/bin/env python

import commands
import os
import glob
import logging
import time
import re 

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.DEBUG )

#challenge_data_path = "/data/celpp/2015/dataset.week.47/stage.2.blastnfilter"
#pdb_protein_location = "/data/pdb.extracted"
#full_copy_location = ""
#target_copy_location = ""

# Name for python binary
PYTHON_BINARY_NAME = 'python'

rdkit_smiles_to_3d_sdf_text = '''
import rdkit.Chem
import rdkit.Chem.AllChem
import sys

if not(len(sys.argv)) == 3:
    print "python smiles2Mol.py inputSmiles outputSdf"
    sys.exit()

smiles = open(sys.argv[1]).read().strip()
mol = rdkit.Chem.MolFromSmiles(smiles)
molH = rdkit.Chem.AddHs(mol)
rdkit.Chem.AllChem.EmbedMolecule(molH)
rdkit.Chem.AllChem.UFFOptimizeMolecule(molH)

w = rdkit.Chem.SDWriter(sys.argv[2])
w.write(molH)
w.close()
'''


chimera_prep_text = '''
import chimera
import sys
opened = chimera.openModels.open(sys.argv[1])
mol = opened[0]

import DockPrep

DockPrep.prep([mol])
from WriteMol2 import writeMol2
with open(sys.argv[2],'wb') as of:
    writeMol2([mol], of)
'''

def extract_info_from_s2(stage_2_out):
    info_dic = {}
    f = open(stage_2_out,"r")
    lines = f.readlines()
    f.close()
    for line in lines:
        #print line.split(",")[0]
        info_title = line.split(",")[0]
        info_value_1 = line.split(",")[1] 
        if line.split(",")[0] == "inchi":
            info_title = line.split(", InChI=")[0] 
            info_value_1 = line.split(", InChI=")[1].split("\n")[0]
        else:
            try:
                info_value_2 = line.split(",")[2].split("\n")[0].strip()
                info_title = line.split(",")[0]
                info_value_1 = line.split(",")[1].strip() 
            except:
                info_value_2 = False
                info_title = line.split(",")[0]
                info_value_1 = line.split(",")[1].split("\n")[0].strip()
                pass
        #if not 
        if not info_value_2:
            #just use the first item if there are multiple entries (use the highest resolution one)
            if info_title not in info_dic:
                info_dic[info_title] = info_value_1
        else:
            if info_title not in info_dic:
                info_dic[info_title] = (info_value_1,info_value_2) 
    return info_dic

#copy all the txt files from the output stage 2 location and create folder for each of this named by the query entry
#check if it finished later needed 

def ligand_prepare(ligand_smile, out_lig_file, rdkit_python):
#    commands.getoutput("$SCHRODINGER/ligprep -WAIT -i 0 -nt -s 1 -g -ismi %s -omae %s"%(ligand_smile, out_lig_file) ) 
#    return os.path.isfile(out_lig_file)
    if os.path.isfile(out_lig_file):
        logging.info('Ligand file %s is already prepared. Skipping.' %(out_lig_file))
        return True

    # Prepare a 3D version of the ligand using babel
    unprep_lig_file_1 = ligand_smile.replace('.smi','_unprep_step1.sdf')
    with open('rdkit_smiles_to_3d_sdf.py','wb') as of:
        of.write(rdkit_smiles_to_3d_sdf_text)

    rdkitpythonpath = os.path.join(rdkit_python, PYTHON_BINARY_NAME)
    commands.getoutput(rdkitpythonpath + ' rdkit_smiles_to_3d_sdf.py ' +
                       ligand_smile + ' ' + unprep_lig_file_1 +
                       ' > rdkit_smiles_to_3d_sdf_out 2>&1')
    unprep_lig_file_2 = ligand_smile.replace('.smi','_unprep_step2.mol2')
    commands.getoutput('babel -isdf %s -omol2 %s' %(unprep_lig_file_1, unprep_lig_file_2))
    #unprep_lig_file = ligand_smile.replace('.smi','_unprep.mol2')
    #commands.getoutput('babel -ismi %s -omol2 %s --gen3D' %(ligand_smile,unprep_lig_file))
    with open('chimeraPrep.py','wb') as of:
        of.write(chimera_prep_text) 
    commands.getoutput('chimera --nogui --script "chimeraPrep.py %s %s" >& chimeraLigPrep.out' %(unprep_lig_file_2, out_lig_file))
    #time.sleep(sleep_time) 
    return os.path.isfile(out_lig_file)

def split_complex (part, complex_file, out_split):
    commands.getoutput("$SCHRODINGER/run split_structure.py -many_files -m %s %s %s"%(part, complex_file, out_split))
    out_receptor = os.path.splitext(out_split)[0] + "_receptor1"+ os.path.splitext(out_split)[1]
    if os.path.isfile(out_receptor):
        return out_receptor 
    else:
        return False

def prepare_protein (protein_file, prepared_protein, sleep_time = 300):
    #here the prepared protein should has ".mol2" as extension
    #commands.getoutput("$SCHRODINGER/utilities/prepwizard %s %s"%(protein_file, prepared_protein))
    with open('chimeraPrep.py','wb') as of:
        of.write(chimera_prep_text) 
    commands.getoutput('chimera --nogui --script "chimeraPrep.py %s %s" >& chimeraProtPrep.out' %(protein_file, prepared_protein))
    #time.sleep(sleep_time) 
    if os.path.isfile(prepared_protein):
        return True
    else:
        return False


def main_proteinprep (challenge_data_path, pdb_protein_path, working_folder, rdkit_python ):
    os.chdir(working_folder)
    current_dir_layer_1 = os.getcwd()

    ## Get all potential target directories and candidates within
    valid_candidates = {}
    #target_ligands = {}
    pot_target_dirs = list(os.walk(challenge_data_path))[0][1]
    #target_ids = []
    # Ensure that the directories are valid
    for pot_target_id in pot_target_dirs:
        os.chdir(current_dir_layer_1)
        # Does it look like a pdb id?
        if len(pot_target_id) != 4:
            logging.info('Filtering potential target directories: %s is not 4 characters long. Skipping' %(pot_target_id))
            continue
        commands.getoutput('mkdir %s' %(pot_target_id))
        #target_ids.append(pot_target_dir)
        valid_candidates[pot_target_id] = []
        target_dir_path = os.path.join(challenge_data_path, pot_target_id)

        # Pull in the ligand inchi
        lig_smiles_files = glob.glob('%s/lig_*.smi' %(target_dir_path))
        if len(lig_smiles_files) != 1:
            logging.info('Unable to find unambiguous ligand smiles for %s - glob returned %r' %(pot_target_id, lig_smiles_files))
            continue
        lig_smiles_file = lig_smiles_files[0]
        local_smiles_file = os.path.basename(lig_smiles_file)
        dest_smiles_file = os.path.join(pot_target_id, local_smiles_file)
        commands.getoutput('cp %s %s' %(lig_smiles_file, dest_smiles_file))
        center_file = os.path.join(target_dir_path,'center.txt')
        commands.getoutput('cp %s %s' %(center_file, pot_target_id))
        ## Get the ligand center of mass for the "LMCSS" candidate (all of the other candidates have been aligned to this one)
        LMCSS_ligand_filenames = glob.glob('%s/LMCSS-*-lig.pdb'%(target_dir_path))
        if len(LMCSS_ligand_filenames) != 1:
            logging.info("Failed to find LMCSS structure's ligand file. There should be one match but I found %r" %(LMCSS_ligand_filenames))
        LMCSS_ligand_filename = LMCSS_ligand_filenames[0]
        
        # Copy in each valid candidate
        for candidate_file in glob.glob('%s/*-%s_*.pdb' %(target_dir_path, pot_target_id)):
            # The LMCSS ligand will be in a pdb file called something like celpp_week19_2016/1fcz/LMCSS-1fcz_1fcz-156-lig.pdb
            # We want to make sure we don't treat this like a receptor
            if 'lig.pdb' in candidate_file:
                continue
            commands.getoutput('cp %s %s' %(candidate_file, pot_target_id))
            candidate_local_file = os.path.basename(candidate_file)
            valid_candidates[pot_target_id].append((local_smiles_file, candidate_local_file))
                
    for target_id in valid_candidates.keys():
        os.chdir(target_id)

        ######################
        #step 6, prepare all proteins
        ######################
        for smiles_filename, candidate_filename in valid_candidates[target_id]:
            ## Parse the candidate name 
            ## Get the method type, target, and candidate info from the filename
            # for example, this will parse 'hiResApo-5hib_2eb2_docked.mol' into [('hiResApo', '5hib', '2eb2')]
            
            parsed_name = re.findall('([a-zA-Z0-9]+)-([a-zA-Z0-9]+)_([a-zA-Z0-9]+)-?([a-zA-Z0-9]*).pdb', candidate_filename)
            if len(parsed_name) != 1:
                logging.info('Failed to parse docked structure name "%s". Parsing yielded %r' %(candidate_filename, parsed_name))
                continue
            candidate_structure_type = parsed_name[0][0]
            candidate_structure_target = parsed_name[0][1]
            candidate_structure_candidate = parsed_name[0][2]
            candidate_structure_ligand = parsed_name[0][2]
            
            # Prepare the ligand
            if not ligand_prepare(smiles_filename, smiles_filename.replace('.smi','_prepared.mol2'), rdkit_python):
                logging.info("Unable to prepare the ligand for this query protein:%s"%target_id)
                #os.chdir(current_dir_layer_1)
                continue 

            # Split the complex 
            #candidate_prefix = candidate_filename.replace('.pdb','')
            candidate_prefix = '%s-%s_%s' %(candidate_structure_type,
                                            candidate_structure_target,
                                            candidate_structure_candidate)
            
            out_split = candidate_prefix+ "_split.pdb"
            out_receptor = split_complex("pdb", candidate_filename, out_split)
            
            
            if not out_receptor:
                logging.info("Unable to split this protein:%s"%(candidate_filename))
                #os.chdir(current_dir_layer_1)
                continue

            logging.info("Successfully split this protein:%s, go to preparation step"%(candidate_filename))
            prepared_protein_mol2 = candidate_prefix + "_prepared.mol2"
            #pass the wizard sleep time here
            preparation_result = prepare_protein(out_receptor,prepared_protein_mol2, 180 )
            if not preparation_result:
                logging.info("Unable to prepare this protein:%s"%(out_split))
                #os.chdir(current_dir_layer_1)
                continue
            #convert into pdb format
            out_prepare_pdb = candidate_prefix + "_prepared.pdb"
            commands.getoutput("$SCHRODINGER/utilities/pdbconvert -imae %s -opdb %s"%(prepared_protein_mol2, out_prepare_pdb))
            logging.info("Successfully prepared this protein:%s"%(out_prepare_pdb))
        os.chdir(current_dir_layer_1)
                        
                

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-p", "--pdbdb", metavar = "PATH", help = "PDB DATABANK which we will dock into")
    parser.add_argument("-c", "--candidatedir", metavar="PATH", help = "PATH where we could find the stage 2 output")
    parser.add_argument("-o", "--outdir", metavar = "PATH", help = "PATH where we run stage 3")
    parser.add_argument("-r", "--rdkitpython", metavar="PATH",
                        help="Path for python build with new "
                             "version of rdkit.",
                        default="")
    #parser.add_option("-s", "--sleep", metavar = "VALUE", help = "Sleep time for protein prep")
    #parser.add_option("-u", "--update", default = False, action = "store_true", help = "update the protein generation and docking step")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.DEBUG )
    opt = parser.parse_args()
    pdb_location = opt.pdbdb
    challenge_data_path = opt.candidatedir
    result_path = opt.outdir
    rdkit_python = opt.rdkitpython
    #sleep_time = opt.sleep
    #running under this dir
    running_dir = os.getcwd()
    main_proteinprep(challenge_data_path, pdb_location, result_path, rdkit_python)
    #move the final log file to the result dir
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s"%(log_file_path, result_path))

