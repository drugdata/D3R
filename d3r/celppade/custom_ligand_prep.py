#!/usr/bin/env python

__author__ = 'j5wagner'

import commands
import os
import glob
import logging
import time
import re 
import shutil
from d3r.utilities.challenge_data import ChallengeData

logger = logging.getLogger(__name__)



class LigandPrep(object):
    """Abstract class defining methods for a custom protein and ligand preparation solution
    for CELPP
    """



    # This prep script will be required to output files with the appropriate suffixes
    OUTPUT_LIG_SUFFIX = '.smi'
    


    
    def ligand_scientific_prep(self, lig_smi_file, out_lig_file, info_dic={}):
       """Does not do any scientific preparation - Passes ligand smiles file forward without any processing 
       """
       
       shutil.copyfile(lig_smi_file, out_lig_file)
       return True
    


    def run_scientific_ligand_prep(self, challenge_data_path, pdb_protein_path, working_folder):
        abs_challenge_data_path = os.path.abspath(challenge_data_path)
        chal_data_obj = ChallengeData(abs_challenge_data_path)
        if not(chal_data_obj.is_valid_for_celpp):
            logging.info('%s is not a valid CELPP challenge data directory. Unable to run ligand prep.')
            return False
        week_chal_data_dict = chal_data_obj.get_targets()
        week_name = week_chal_data_dict.keys()[0]
        abs_week_path = os.path.join(abs_challenge_data_path, week_name)
        pot_target_dirs = week_chal_data_dict[week_name]
        os.chdir(working_folder)
        current_dir_layer_1 = os.getcwd()

        ## Get all potential target directories and candidates within
        valid_targets = {}

        # Ensure that the directories are valid
        for pot_target_dir in pot_target_dirs:
            os.chdir(current_dir_layer_1)
            pot_target_id = os.path.basename(pot_target_dir.strip('/'))

            # Does it look like a pdb id?
            if len(pot_target_id) != 4:
                logging.info('Filtering potential target directories: %s is not 4 characters long. Skipping' %(pot_target_id))
                continue
            os.mkdir(pot_target_id)
            target_dir_path = os.path.join(abs_week_path, pot_target_id)

            # Pull in the ligand inchi/smiles
            
            lig_smiles_files = glob.glob('%s/lig_*.smi' %(target_dir_path))
            if len(lig_smiles_files) != 1:
                logging.info('Unable to find unambiguous ligand smiles for %s - glob returned %r' %(pot_target_id, lig_smiles_files))
                continue
            lig_smiles_file = lig_smiles_files[0]
            local_smiles_file = os.path.basename(lig_smiles_file)
            dest_smiles_file = os.path.join(pot_target_id, local_smiles_file)
            shutil.copyfile(lig_smiles_file, dest_smiles_file)
            
            valid_targets[pot_target_id] = local_smiles_file


            
        for target_id in valid_targets.keys():
            os.chdir(target_id)

            smiles_filename = valid_targets[target_id]

            # Prepare the ligand
            lig_prefix = smiles_filename.replace('.smi','')
            prepared_lig_file = '%s_prepared%s' %(lig_prefix, LigandPrep.OUTPUT_LIG_SUFFIX)
            lig_prep_result = self.ligand_scientific_prep(smiles_filename, prepared_lig_file)
            if lig_prep_result == False:
                logging.info("Unable to prepare the ligand for this target protein: %s. Skipping" %(target_id))
                continue 

            if not(os.path.exists(prepared_lig_file)):
                logging.info('Expected output file %s does not exist. Assuming that ligand prep failed. Skipping target %s' %(prepared_lig_file, target_id))
                continue
            if os.path.getsize(prepared_lig_file)==0:
                logging.info('Expected output file %s has size 0. Assuming that ligand prep failed. Skipping candidate %s' %(prepared_lig_file, target_id))
                continue

            logging.info("Successfully prepared ligand %s for target %s"%(lig_prefix, target_id))


            os.chdir(current_dir_layer_1)

