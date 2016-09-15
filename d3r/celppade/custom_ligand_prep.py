#!/usr/bin/env python

__author__ = 'j5wagner'

import commands
import os
import glob
import logging
import time
import re 
import shutil

logger = logging.getLogger(__name__)



class LigandPrep(object):
    """Abstract class defining methods for a custom protein and ligand preparation solution
    for CELPP
    """



    # This prep script will be required to output files with the appropriate suffixes
    '''
    OUTPUT_PROTEIN_SUFFIX = '.pdb'
    '''
    OUTPUT_LIG_SUFFIX = '.smi'
    


    
    def ligand_scientific_prep(self, lig_smi_file, out_lig_file, info_dic={}):
       """Does not do any scientific preparation - Passes ligand smiles file forward without any processing 
        """
        
        shutil.copyfile(lig_smi_file, out_lig_file)
        return True
    
    '''
    def prepare_protein(self, protein_file, prepared_protein_file, info_dic={}):
        """Does not do any scientific preparation - Passes protein forward without any processing 
        """
 
        commands.getoutput('cp %s %s' %(protein_file, prepared_protein_file))
        if not(os.path.isfile(prepared_protein_file)):
            logging.info('Prepared protein file %s does not exist. Assuming that prep failed.' %(prepared_protein_file))
            return False
        else:
            return True
            


    def split_complex(self, complex_pdb):
        # This function takes a pdb name of (potentially) a protein
        # complex, and should return only the atoms that should be
        # seen during docking

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

    '''




    def run_scientific_ligand_prep(self, challenge_data_path, pdb_protein_path, working_folder):
        abs_challenge_data_path = os.path.abspath(challenge_data_path)
        os.chdir(working_folder)
        current_dir_layer_1 = os.getcwd()

        ## Get all potential target directories and candidates within
        valid_targets = {}
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
            #commands.getoutput('mkdir %s' %(pot_target_id))
            os.mkdir(pot_target_id)
            #target_ids.append(pot_target_dir)
            '''
            valid_targets[pot_target_id] = []
            '''
            target_dir_path = os.path.join(abs_challenge_data_path, pot_target_id)

            # Pull in the ligand inchi/smiles
            
            lig_smiles_files = glob.glob('%s/lig_*.smi' %(target_dir_path))
            if len(lig_smiles_files) != 1:
                logging.info('Unable to find unambiguous ligand smiles for %s - glob returned %r' %(pot_target_id, lig_smiles_files))
                continue
            lig_smiles_file = lig_smiles_files[0]
            local_smiles_file = os.path.basename(lig_smiles_file)
            dest_smiles_file = os.path.join(pot_target_id, local_smiles_file)
            shutil.copyfile(lig_smiles_file, dest_smiles_file))
            
            '''
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
                '''
            valid_targets[pot_target_id] = local_smiles_file
            
        for target_id in valid_targets.keys():
            os.chdir(target_id)

            smiles_filename = valid_targets[target_id]
            '''
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
            '''
            # Prepare the ligand
            lig_prefix = smiles_filename.replace('.smi','')
            prepared_lig_file = '%s_prepared%s' %(lig_prefix, Prep.OUTPUT_LIG_SUFFIX)
            lig_prep_result = self.prepare_ligand(smiles_filename, prepared_lig_file)
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

            '''
            # Split the complex 
            candidate_prefix = '%s-%s_%s' %(candidate_structure_type,
                                            candidate_structure_target,
                                            candidate_structure_candidate)

            split_intermediate_prefix = candidate_prefix+ "_split.pdb"
            split_receptor_file = self.split_complex(candidate_filename)


            if not(split_receptor_file):
                logging.info("Unable to split this protein:%s"%(candidate_filename))
                continue

            logging.info("Successfully split this protein:%s, going to preparation step"%(candidate_filename))
            prepared_protein_file = "%s_prepared%s" %(candidate_prefix, Prep.OUTPUT_PROTEIN_SUFFIX)

            preparation_result = self.prepare_protein(split_receptor_file, prepared_protein_file)
            if not preparation_result:
                logging.info("Unable to prepare this protein:%s"%(split_receptor_file))
                continue
            
            #convert into pdb format
            logging.info("Successfully prepared this protein:%s"%(prepared_protein_file))
            '''
            
            os.chdir(current_dir_layer_1)

