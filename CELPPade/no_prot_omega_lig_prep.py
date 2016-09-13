#!/usr/bin/env python

__author__ = 'j5wagner'

import commands
import os
import glob
import logging
import time
import re 

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.DEBUG )







################################################################
########## CONTESTANT MODIFICATIONS SHOULD BEGIN HERE ##########
################################################################



# This prep script will be required to output files with the appropriate suffixes
output_protein_suffix = '.pdb'
output_lig_suffix = '.sdf'



def ligand_prepare(lig_smi_file, out_lig_file, info_dic={}):
    
    lig_prefix = os.path.basename(lig_smi_file).replace('.smi','')
    
    # Perform conformer generation using omega
    omega_stdout_file = lig_prefix + '_omega_confgen_stdout'
    omega_stderr_file = lig_prefix + '_omega_confgen_stderr'
    omega_cmd = 'omega2 -flipper true -in ' + lig_smi_file + ' -out ' + out_lig_file + ' 2> ' + omega_stderr_file + ' 1> ' + omega_stdout_file
    logging.info('Running omega command: ' + omega_cmd)
    commands.getoutput(omega_cmd)
    

    if not(os.path.isfile(out_lig_file)):
        logging.info('Prepared ligand file %s does not exist. Assuming that prep failed.' %(out_lig_file))
        return False
    else:
        return True


def prepare_protein (protein_file, prepared_protein_file, info_dic={}):
    # We don't do any scientific prep here. Openeye docking can handle unprepared proteins
    commands.getoutput('cp %s %s' %(protein_file, prepared_protein_file))
    if not(os.path.isfile(prepared_protein_file)):
        logging.info('Prepared protein file %s does not exist. Assuming that prep failed.' %(prepared_protein_file))
        return False
    else:
        return True



################################################################
########### CONTESTANT MODIFICATIONS SHOULD END HERE ###########
################################################################




def split_complex (complex_pdb):
    # This function takes a pdb name of (potentially) a protein
    # complex, and should return only the relevant protein part for
    # docking
    
    # Poor man's receptor splitting - Remove all HETATM and
    # CONECT lines
    data = open(complex_pdb).readlines()
    output_file = complex_pdb.replace('.pdb','_split.pdb')
    with open(output_file,'wb') as of:
        for line in data:
            if (line[:6] == 'HETATM') or (line[:6] == 'CONECT'):
                continue
            else:
                of.write(line)

    # Ensure that the new file exists
    if not(os.path.isfile(output_file)):
        logging.info('Split pdb file %s not created' %(output_file))
        return False

    # Ensure that the new file contains any data
    if os.path.getsize(output_file) == 0:
        logging.info('Splitting pdb file %s yielded empty file %s' %(complex_file, output_file))
        return False

    return output_file




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



def main_proteinprep (challenge_data_path, pdb_protein_path, working_folder):
    abs_challenge_data_path = os.path.abspath(challenge_data_path)
    os.chdir(working_folder)
    current_dir_layer_1 = os.getcwd()
 
    ## Get all potential target directories and candidates within
    valid_candidates = {}
    #target_ligands = {}
    pot_target_dirs = glob.glob('%s/????' %(abs_challenge_data_path))
    #target_ids = []
    # Ensure that the directories are valid
    for pot_target_dir in pot_target_dirs:
        os.chdir(current_dir_layer_1)
        pot_target_id = os.path.basename(pot_target_dir.strip('/'))
        # Does it look like a pdb id?
        if len(pot_target_id) != 4:
            logging.info('Filtering potential target directories: %s is not 4 characters long. Skipping' %(pot_target_id))
            continue
        commands.getoutput('mkdir %s' %(pot_target_id))
        #target_ids.append(pot_target_dir)
        valid_candidates[pot_target_id] = []
        target_dir_path = os.path.join(abs_challenge_data_path, pot_target_id)

        # Pull in the ligand inchi/smiles
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
            lig_prefix = smiles_filename.replace('.smi','')
            prepared_lig_file = '%s_prepared%s' %(lig_prefix, output_lig_suffix)
            lig_prep_result = ligand_prepare(smiles_filename, prepared_lig_file)
            if lig_prep_result == False:
                logging.info("Unable to prepare the ligand for this target protein: %s. Skipping" %(target_id))
                continue 

            # Split the complex 
            candidate_prefix = '%s-%s_%s' %(candidate_structure_type,
                                            candidate_structure_target,
                                            candidate_structure_candidate)
            
            split_intermediate_prefix = candidate_prefix+ "_split.pdb"
            split_receptor_file = split_complex(candidate_filename)
            
            
            if not(split_receptor_file):
                logging.info("Unable to split this protein:%s"%(candidate_filename))
                continue
            
            logging.info("Successfully split this protein:%s, going to preparation step"%(candidate_filename))
            prepared_protein_file = "%s_prepared%s" %(candidate_prefix, output_protein_suffix)

            preparation_result = prepare_protein(split_receptor_file, prepared_protein_file)
            if not preparation_result:
                logging.info("Unable to prepare this protein:%s"%(split_receptor_file))
                continue
            #convert into pdb format
            logging.info("Successfully prepared this protein:%s"%(prepared_protein_file))


        os.chdir(current_dir_layer_1)
                        
                

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-p", "--pdbdb", metavar = "PATH", help = "PDB DATABANK which we will dock into")
    parser.add_argument("-c", "--challengedata", metavar="PATH", help = "PATH to the unpacked challenge data package")
    parser.add_argument("-o", "--prepdir", metavar = "PATH", help = "PATH to the output directory")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level = logging.INFO )
    opt = parser.parse_args()
    pdb_location = opt.pdbdb
    challenge_data_path = opt.challengedata
    prep_result_path = opt.prepdir

    #running under this dir
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(prep_result_path), 'final.log')

    main_proteinprep(challenge_data_path, pdb_location, prep_result_path)

    #move the final log file to the result dir
    commands.getoutput("mv %s %s"%(log_file_path, prep_result_path))

