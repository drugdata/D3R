#!/usr/bin/env python

__author__ = 'j5wagner'

import commands
import sys
import os
import glob
import logging
import time
import re 
import shutil 
from d3r.utilities.challenge_data import ChallengeData
from d3r.utilities.readers import ReadText


logger = logging.getLogger(__name__)



class ProteinPrep(object):
    """Abstract class defining methods for a custom protein and ligand preparation solution
    for CELPP
    """



    # This prep script will be required to output files with the appropriate suffixes
    OUTPUT_PROTEIN_SUFFIX = '.pdb'

    def receptor_scientific_prep(self, protein_file, prepared_protein_file, targ_info_dict={}):
        """Does not do any scientific preparation - Passes protein forward without any processing 
        """
        
        # Poor man's receptor splitting - Remove all HETATM and
        # CONECT lines
        #data = open(protein_file).readlines()
        #with open(prepared_protein_file,'wb') as of:
        #    for line in data:
        #        if (line[:6] == 'HETATM') or (line[:6] == 'CONECT'):
        #            continue
        #        else:
        #            of.write(line)
        
        shutil.copyfile(protein_file, prepared_protein_file)
        return True





    def run_scientific_protein_prep(self, challenge_data_path, pdb_protein_path, working_folder):
        abs_challenge_data_path = os.path.abspath(challenge_data_path)
        chal_data_obj = ChallengeData(abs_challenge_data_path)
        if not(chal_data_obj.is_valid_for_celpp()):
            logging.info('%s is not a valid CELPP challenge data directory. Unable to run protein prep.' %(abs_challenge_data_path))
            return False
        week_chal_data_dict = chal_data_obj.get_targets()
        week_name = week_chal_data_dict.keys()[0]
        abs_week_path = os.path.join(abs_challenge_data_path, week_name)
        pot_target_dirs = week_chal_data_dict[week_name]
        os.chdir(working_folder)
        current_dir_layer_1 = os.getcwd()

        ## Get all potential target directories and candidates within
        valid_candidates = {}


        # Ensure that the chalengedata targets are valid and copy in files
        for pot_target_dir in pot_target_dirs:
            os.chdir(current_dir_layer_1)
            pot_target_id = os.path.basename(pot_target_dir.strip('/'))
            # Does it look like a pdb id?
            if len(pot_target_id) != 4:
                logging.info('Filtering potential target directories: %s is not 4 characters long. Skipping' %(pot_target_id))
                continue

            os.mkdir(pot_target_id)

            valid_candidates[pot_target_id] = []
            target_dir_path = os.path.join(abs_week_path, pot_target_id)

            # Copy in <targ id>.txt file
            targ_info_basename = pot_target_id+'.txt'
            targ_info_file = os.path.join(target_dir_path, targ_info_basename)
            targ_info_dest = os.path.join(pot_target_id, targ_info_basename)
            shutil.copyfile(targ_info_file, targ_info_dest)
            
            # Copy in center.txt file
            center_file = os.path.join(target_dir_path,'center.txt')
            #center_file_basename = os.path.basename(center_file)
            center_file_dest = os.path.join(pot_target_id, 'center.txt')
            shutil.copyfile(center_file, center_file_dest)


            # Copy in each valid candidate
            for candidate_file in glob.glob('%s/*-%s_*.pdb' %(target_dir_path, pot_target_id)):
                # The LMCSS ligand will be in a pdb file called something like celpp_week19_2016/1fcz/LMCSS-1fcz_1fcz-156-lig.pdb
                # We want to make sure we don't treat this like a receptor
                if 'lig.pdb' in candidate_file:
                    continue
                candidate_file_basename = os.path.basename(candidate_file)
                candidate_file_dest = os.path.join(pot_target_id,candidate_file_basename)
                shutil.copyfile(candidate_file, candidate_file_dest)
                candidate_local_file = os.path.basename(candidate_file)
                valid_candidates[pot_target_id].append(candidate_local_file)

        for target_id in valid_candidates.keys():
            os.chdir(current_dir_layer_1)
            os.chdir(target_id)
            current_dir_layer_2 = os.getcwd()          
            ReadText_obj = ReadText()
            targ_info_dict = ReadText_obj.parse_txt(target_id + '.txt')
            

            for candidate_filename in valid_candidates[target_id]:
                os.chdir(current_dir_layer_2)
                ## Parse the candidate name 
                ## Get the method type, target, and candidate info from the filename
                # for example, this will parse 'hiResApo-5hib_2eb2.pdb' into [('hiResApo', '5hib', '2eb2')]

                parsed_name = re.findall('([a-zA-Z0-9]+)-([a-zA-Z0-9]+)_([a-zA-Z0-9]+)-?([a-zA-Z0-9]*).pdb', candidate_filename)
                if len(parsed_name) != 1:
                    logging.info('Failed to parse docked structure name "%s". Parsing yielded %r' %(candidate_filename, parsed_name))
                    continue
                candidate_structure_type = parsed_name[0][0]
                candidate_structure_target = parsed_name[0][1]
                candidate_structure_candidate = parsed_name[0][2]
                candidate_structure_ligand = parsed_name[0][2]

                
                candidate_prefix = '%s-%s_%s' %(candidate_structure_type,
                                                candidate_structure_target,
                                                candidate_structure_candidate)          
                # Make candidate prep directory
                os.mkdir(candidate_prefix)
                # Copy in raw candidate file
                candidate_copy_origin = candidate_filename
                candidate_copy_dest = os.path.join(candidate_prefix, candidate_filename)
                shutil.copyfile(candidate_copy_origin, candidate_copy_dest)
                # Copy in center file
                center_copy_origin = 'center.txt'
                center_copy_dest = os.path.join(candidate_prefix, 'center.txt')
                shutil.copyfile(center_copy_origin, center_copy_dest)
                # Move into candidate prep directory
                os.chdir(candidate_prefix)

                # Run prep
                prepared_protein_file = "%s_prepared%s" %(candidate_prefix, ProteinPrep.OUTPUT_PROTEIN_SUFFIX)

                try:
                    preparation_result = self.receptor_scientific_prep(candidate_filename, 
                                                                       prepared_protein_file,
                                                                       targ_info_dict=targ_info_dict)
                except:
                    logging.info(sys.exc_info())
                    logging.info('try/except statement caught error in scientific protein prep. Skipping candidate %s' %(candidate_prefix))
                    continue

                if preparation_result == False:
                    logging.info("Unable to prepare this protein:%s"%(candidate_filename))
                    continue
                if not(os.path.exists(prepared_protein_file)):
                    logging.info('Expected output file %s does not exist after protein preparation. Assuming that protein prep failed. Skipping candidate %s' %(prepared_protein_file, candidate_prefix))
                    continue
                if os.path.getsize(prepared_protein_file)==0:
                    logging.info('Expected output file %s has size 0. Assuming that protein prep failed. Skipping candidate %s' %(prepared_protein_file, candidate_prefix))
                    continue

                prepared_receptor_origin = prepared_protein_file
                prepared_receptor_dest = os.path.join(current_dir_layer_2, prepared_protein_file)
                shutil.copyfile(prepared_receptor_origin, prepared_receptor_dest)          
                
                logging.info("Successfully prepared this protein:%s"%(prepared_protein_file))


            os.chdir(current_dir_layer_1)

