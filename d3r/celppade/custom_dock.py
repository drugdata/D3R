#!/usr/bin/env python

__author__ = 'j5wagner'

import os
import sys
import shutil
import re
import glob
import logging
from d3r.utilities.readers import ReadText

logger = logging.getLogger(__name__)

class Dock(object):
    """Abstract class defining methods for a custom docking solution
    for CELPP
    """

    ## To ensure that the workflow copies over the correct scientific prep
    ## outputs, list the expected filetypes here. This script will look
    ## for files in the target prep directory that follow the naming
    ## convention "<cand category>-<target pdbid>_<cand
    ## pdbid>_prepared.<sci_prepped_prot_suffix>", for example
    ## "LMCSS-5b4t_3w8e_prepared.pdb". It will then pass the appropriate
    ## filenames and prefixes to the technical_prep() and dock()
    ## functions.
    SCI_PREPPED_LIG_SUFFIX = '_prepared.sdf'
    SCI_PREPPED_PROT_SUFFIX = '_prepared.pdb'

    def ligand_technical_prep(self, 
                              sci_prepped_lig, 
                              targ_info_dict={}):
        """
        'Technical preparation' is the step immediate preceding
        docking. During this step, you may perform any file
        conversions or processing that are specific to your docking
        program. Implementation of this function is optional.
        :param sci_prepped_lig: Scientifically prepared ligand file
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: A list of result files to be copied into the
        subsequent docking folder. The base implementation merely
        returns the input string in a list (ie. [sci_prepped_lig]) 
        """
        return [sci_prepped_lig]
    
    def receptor_technical_prep(self, 
                                sci_prepped_receptor, 
                                pocket_center, 
                                targ_info_dict={}):
        """
        'Technical preparation' is the step immediately preceding
        docking. During this step, you may perform any file
        conversions or processing that are specific to your docking
        program. Implementation of this function is optional.
        :param sci_prepped_receptor: Scientifically prepared receptor file
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: A list of result files to be copied into the
        subsequent docking folder. This implementation merely
        returns the input string in a list (ie [sci_prepped_receptor])
        """
        return [sci_prepped_receptor]

    def dock(self, 
             tech_prepped_lig_list, 
             tech_prepped_receptor_list, 
             output_receptor_pdb, 
             output_lig_mol,
             targ_info_dict={}):
        """
        This function is the only one which the contestant MUST
        implement.  The dock() step runs the actual docking
        algorithm. Its first two arguments are the return values from
        the technical preparation functions for the ligand and
        receptor. These arguments are lists of file names (strings),
        which can be assumed to be in the current directory. 
        If prepare_ligand() and ligand_technical_prep() are not
        implemented by the contestant, tech_prepped_lig_list will
        contain a single string which names a SMILES file in the
        current directory.
        If prepare_protein() and receptor_technical_prep() are not
        implemented by the contestant, tech_prepped_receptor_list will
        contain a single string which names a PDB file in the current
        directory.
        The outputs from this step must be two files - a pdb with the
        filename specified in the output_receptor_pdb argument, and a
        mol with the filename specified in the output_ligand_mol
        argument.
        :param tech_prepped_lig_list: The list of file names resturned by ligand_technical_prep. These have been copied into the current directory.
        :param tech_prepped_receptor_list: The list of file names resturned by receptor_technical_prep. These have been copied into the current directory.
        :param output_receptor_pdb: The final receptor (after docking) must be converted to pdb format and have exactly this file name.
        :param output_lig mol: The final ligand (after docking) must be converted to mol format and have exactly this file name.
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: True if docking is successful, False otherwise. Unless overwritten, this implementation always returns False
        """
        return False


    def get_pocket_center(self, target_prep_dir):
        pocket_center_file = os.path.join(target_prep_dir,'center.txt')
        try:
            pocket_center = open(pocket_center_file,"r").readlines()[0]
        except:
            logging.info("Unable to find the center file for this case %s" %(target_prep_dir))
            return False
        pocket_center = pocket_center.split(',')
        if not len(pocket_center) == 3:
            logging.info("Pocket center does not parse into 3 coordinates for case %s. Pocket center is %r" %(target_prep_dir, pocket_center))
            return False

        try:
            pocket_center = [float(i.strip()) for i in pocket_center]
        except:
            logging.info('Error converting pocket center coordinates to float for case %s. Pocket center is %r.' %(target_prep_dir, pocket_center))
            return False
        return pocket_center

    def get_sci_prepped_lig(self, target_prep_dir, sci_prepped_lig_suffix):
        sci_prepped_lig_files = glob.glob('%s/lig_*%s' %(target_prep_dir, sci_prepped_lig_suffix))
        # Ensure that there is one unambiguous ligand for docking
        if len(sci_prepped_lig_files) == 0:
            logging.info('No ligand files found for target prep dir %s.' %(target_prep_dir))
            return False
        if len(sci_prepped_lig_files) > 1:
            logging.info('Multiple ligand files found for target prep dir %s (found ligands %r). The workflow currently should only be sending one ligand. Skipping.' %(target_prep_dir, sci_prepped_lig_files))
            return False
        sci_prepped_lig_file = sci_prepped_lig_files[0]
        return sci_prepped_lig_file

    def parse_lig_filename(self, sci_prepped_lig_file):
        try:
            sci_prepped_lig_basename = os.path.basename(sci_prepped_lig_file)
            lig_re_pattern = 'lig_([a-zA-Z0-9]{3})%s' %(self.SCI_PREPPED_LIG_SUFFIX)
            lig_re_results = re.findall(lig_re_pattern, sci_prepped_lig_basename)
        except:
            logging.info('Error parsing ligand filename %s. Skipping.' %(sci_prepped_lig_file))
            return False
        if len(lig_re_results) != 1:
            logging.info('Unable to parse prepared ligand %s. Regular expression matching yielded %r' %(sci_prepped_lig_file, lig_re_results))
            return False
        lig_prefix = lig_re_results[0]
        return lig_prefix

    def parse_cand_name(self, cand_name):
        try:
            cand_re_pattern = '([a-zA-Z0-9]+)-([a-zA-Z0-9]{4})_([a-zA-Z0-9]{4})%s' %(self.SCI_PREPPED_PROT_SUFFIX)
            cand_re_results = re.findall(cand_re_pattern, cand_name)
        except:
            logging.info('Error parsing prepared protein %s. Skipping' %(cand_name))
            return False
        if len(cand_re_results) != 1:
            logging.info('Unable to parse prepared protein %s. Regular expression matching did not yield one match (yielded %r)' %(cand_name, cand_re_results))
            return False
                #raise Exception('Unable to parse prepared protein %s. Regular expression matching yielded %r' %(cand_name, cand_re_results))
        category, target, cand = cand_re_results[0]
        return category, target, cand



    def run_dock(self, prot_sci_prep_dir, lig_sci_prep_dir, dock_dir):
        #os.chdir(prep_result_dir)
        abs_lig_sci_prep_dir = os.path.abspath(lig_sci_prep_dir)
        abs_prot_sci_prep_dir = os.path.abspath(prot_sci_prep_dir)
        abs_dock_dir = os.path.abspath(dock_dir)

        targ_prot_prep_dirs = glob.glob('%s/????' %(abs_prot_sci_prep_dir))
        #targ_lig_prep_dirs = glob.glob('%s/????' %(abs_lig_prep_dir))
        #prepped_prot_targs = [os.path.basename(i.rstrip('/')) for i in targ_prot_prep_dirs]
        #prepped_lig_targs = [os.path.basename(i.rstrip('/')) for i in targ_lig_prep_dirs]
        
        targ_dic = {}
    
        for targ_prot_prep_dir in targ_prot_prep_dirs:
            targ_name = os.path.basename(targ_prot_prep_dir.rstrip('/'))
            targ_lig_prep_dir = os.path.join(abs_lig_sci_prep_dir, targ_name)
            targ_dock_dir = os.path.join(abs_dock_dir,targ_name)
            #os.mkdir(target_dock_dir)
            logging.info("============= Starting to process target:%s =============" %targ_name)
            targ_dic[targ_name] = {}
            # Until we find a good ligand and cand receptor, mark this as invalid
            targ_dic[targ_name]['valid_targ'] = False
        
            # Get the binding pocket center
            pocket_center =  self.get_pocket_center(targ_prot_prep_dir)
            if pocket_center == False:
                logging.info('Failed to find pocket center file in dirctory %s. Skipping target %s.' %(targ_prot_prep_dir, targ_name))
                continue

            # Get the ligand name in this directory
            sci_prepped_lig_file = self.get_sci_prepped_lig(targ_lig_prep_dir,
                                                            Dock.SCI_PREPPED_LIG_SUFFIX)
            if sci_prepped_lig_file == False:
                logging.info('Unable to find single ligand for target dir %s. Skipping target %s.' %(targ_lig_prep_dir, targ_name))
                continue

            # Process ligand file names
            lig_prefix = self.parse_lig_filename(sci_prepped_lig_file)
            if lig_prefix == False:
                logging.info('Unable to parse ligand filename %s. Skipping target %s.' %(sci_prepped_lig_file, targ_name))
                continue

            # Get the cand protein names in this directory
            potential_cand_proteins = glob.glob('%s/*-????_????%s' %(targ_prot_prep_dir, Dock.SCI_PREPPED_PROT_SUFFIX))

            # Process potential cand protein names
            for potential_cand_protein in potential_cand_proteins:
                potential_cand_basename = os.path.basename(potential_cand_protein)
                category, targ_id, cand_id = self.parse_cand_name(potential_cand_basename)
                targ_txt_file = os.path.join(targ_lig_prep_dir,targ_id + '.txt')
                # If this is the first valid cand for this target, make a new dictionary key and mark the target as valid
                if not('valid_cands' in targ_dic[targ_name].keys()):
                    targ_dic[targ_name]['valid_targ'] = True
                    targ_dic[targ_name]['prot_prep_dir'] = targ_prot_prep_dir
                    targ_dic[targ_name]['lig_prep_dir'] = targ_lig_prep_dir
                    targ_dic[targ_name]['pocket_center'] = pocket_center
                    targ_dic[targ_name]['targ_txt_file'] = targ_txt_file
                    targ_dic[targ_name]['lig_file'] = sci_prepped_lig_file
                    targ_dic[targ_name]['lig_prefix'] = lig_prefix
                    targ_dic[targ_name]['valid_cands'] = []
                targ_dic[targ_name]['valid_cands'].append((potential_cand_protein,
                                                          category,
                                                          targ_id,
                                                          cand_id))
            
            
            os.chdir(abs_dock_dir)
    
        # Print out all of the target/cands we will dock to
        logging.info('targ_dic is: ')
        for targ in targ_dic.keys():
            logging.info('%s: %s' %(targ, targ_dic[targ]))

        # Generate a list of targets for docking
        targ_names = targ_dic.keys()
        targ_names = [i for i in targ_names if targ_dic[i]['valid_targ']==True]
        targ_names.sort()

        # Begin populating the target directory
        os.chdir(abs_dock_dir)
        for targ_name in targ_names:
            pocket_center = targ_dic[targ_name]['pocket_center']
            os.chdir(abs_dock_dir)
            # Make the target directory
            os.mkdir(targ_name)
            os.chdir(targ_name)
            abs_targ_dock_dir = os.getcwd()


            # Copy the targ.txt file in 
            copy_dest = '%s/%s' %(abs_targ_dock_dir, os.path.basename(targ_dic[targ_name]['targ_txt_file']))
            shutil.copyfile(targ_dic[targ_name]['targ_txt_file'], copy_dest)
            # Parse the targ.txt file
            ReadText_obj = ReadText()
            targ_info_dict = ReadText_obj.parse_txt(copy_dest)
            targ_info_dict['pocket_center'] = pocket_center

            #### Run CELPPade technical prep

            ### Ligand technical prep

            ## Ligand technical prep setup
            lig_tech_prep_dir = 'lig_%s_tech_prep' %(targ_dic[targ_name]['lig_prefix'])
            os.mkdir(lig_tech_prep_dir)



            # Copy the sci prepped ligand in
            copy_dest = '%s/%s' %(lig_tech_prep_dir, os.path.basename(targ_dic[targ_name]['lig_file']))
            shutil.copyfile(targ_dic[targ_name]['lig_file'], copy_dest)
        
            lig_base_filename = os.path.basename(targ_dic[targ_name]['lig_file'])
            os.chdir(lig_tech_prep_dir)
        
            ## Call user-defined ligand technical prep
            try:
                tech_prepped_lig_file_list = self.ligand_technical_prep(lig_base_filename, targ_info_dict=targ_info_dict)
            except:
                logging.info(sys.exc_info())
                logging.info('try/except statement caught error in function lig_technical_prep. Skipping target %s.' %(targ_name))

                continue
        
            ## Ensure that ligand technical prep was successful
            # Check for function-reported failure
            if tech_prepped_lig_file_list == False:
                logging.info('Technical ligand preparation failed on %s. Skipping target %s.' %(os.path.abspath(lig_base_filename), targ_name))
                continue

            # Ensure that ligand technical prep returns a list of filenames
            if not(type(tech_prepped_lig_file_list) is list):
                logging.info('Technical ligand preparation for %s did not return a list of filenames. Skipping target %s.' %(os.path.abspath(lig_base_filename), targ_name))
                continue

            # Ensure that all files in list really exist
            for filename in tech_prepped_lig_file_list:
                if not(os.path.exists(filename)):
                    logging.info('Technical ligand preparation for %s returned file list %r, but file %s does not exist. Skipping target %s.' %(os.path.abspath(lig_base_filename), tech_prepped_lig_file_list, filename, targ_name))
                    continue
            

            ## Prepare to copy these files for later
            tech_prepped_lig_file_list = [os.path.abspath(i) for i in tech_prepped_lig_file_list]
            targ_dic[targ_name]['tech_prepped_lig_files'] = tech_prepped_lig_file_list

            logging.info('Technical ligand prep successful for %s. All files exist from returned list %r. ' %(os.path.abspath(lig_base_filename), tech_prepped_lig_file_list))

            ### Candidate tech prep
            for cand_file, category, targ_id, cand_id in targ_dic[targ_name]['valid_cands']:
                os.chdir(abs_targ_dock_dir)
            
                ## Candidate tech prep setup
                cand_tech_prep_dir = '%s_%s_tech_prep/' %(category, cand_id)
                os.mkdir(cand_tech_prep_dir)
                copy_dest = '%s/%s' %(cand_tech_prep_dir, os.path.basename(cand_file))
                shutil.copyfile(cand_file, copy_dest)
                prot_base_filename = os.path.basename(cand_file)
                os.chdir(cand_tech_prep_dir)
            
                ## Call user-defined protein technical prep
                try:
                    tech_prepped_prot_file_list = self.receptor_technical_prep(prot_base_filename, pocket_center, targ_info_dict=targ_info_dict)
                except:
                    logging.info(sys.exc_info())
                    logging.info('try/except statement caught error in function receptor_technical_prep.  Skipping candidate %s for target %s.' %(cand_id, targ_name))
                    continue

                ## Ensure that receptor technical prep was successful
                # Check for function-reported failure
                if tech_prepped_prot_file_list == False:
                    logging.info('Technical protein preparation failed on %s. Skipping candidate %s for target %s.' %(os.path.abspath(prot_base_filename), cand_id, targ_name))
                    continue

                # Ensure that receptor technical prep returns a list of filenames
                if not(type(tech_prepped_prot_file_list) is list):
                    logging.info('Technical protein preparation for %s did not return a list of filenames. Skipping candidate %s for target %s.' %(os.path.abspath(prot_base_filename), cand_id, targ_name))
                    continue

                # Ensure that all files in list really exist
                for filename in tech_prepped_prot_file_list:
                    if not(os.path.exists(filename)):
                        logging.info('Technical protein preparation for %s returned file list %r, but file %s does not exist. Skipping candidate %s for target %s.' %(os.path.abspath(lig_base_filename), tech_prepped_lig_file_list, filename, cand_id, targ_name))
                        continue
            
                ## Prepare to copy these files for docking step
                tech_prepped_prot_file_list = [os.path.abspath(i) for i in tech_prepped_prot_file_list]

                logging.info('Protein technical prep successful for %s. All files exist from returned list %r.'%(cand_file, tech_prepped_prot_file_list))
            
        
                #### Run CELPPade docking

                os.chdir(abs_targ_dock_dir)
                cand_dock_dir = '%s_%s_docking' %(category, cand_id)
                os.mkdir(cand_dock_dir)
                os.chdir(cand_dock_dir)

                ## Prepare expected file names
                output_receptor_pdb = '%s-%s_%s_docked.pdb' %(category,
                                                              targ_id,
                                                              cand_id)
                output_lig_mol = '%s-%s_%s_docked.mol' %(category,
                                                         targ_id,
                                                         cand_id)

                ## Copy in tech prepped files
                for filename in tech_prepped_lig_file_list:
                    file_base_name = os.path.basename(filename)
                    shutil.copyfile(filename,file_base_name)
                for filename in tech_prepped_prot_file_list:
                    file_base_name = os.path.basename(filename)
                    shutil.copyfile(filename,file_base_name)
                
                ## Do the actual docking
                try:
                    dock_results = self.dock(tech_prepped_lig_file_list,
                                             tech_prepped_prot_file_list,
                                             output_receptor_pdb,
                                             output_lig_mol,
                                             targ_info_dict=targ_info_dict)
                except:
                    logging.info(sys.exc_info())
                    logging.info('try/except statement caught error in dock() function. Docking was given '
                                 'inputs tech_prepped_lig_file_list=%r tech_prepped_prot_file_list=%r '
                                 'output_receptor_pdb=%r output_lig_mol=%r. Skipping docking to this '
                                 'candidate.'
                                 %(tech_prepped_lig_file_list,
                                   tech_prepped_prot_file_list,
                                   output_receptor_pdb, output_lig_mol))
                    continue

                ## Check for success
                # Check for self-reported failure
                if dock_results == False:
                    logging.info('Docking returned False given inputs: tech_prepped_lig_file_list=%r   '
                                 'tech_prepped_prot_file_list=%r    output_receptor_pdb=%r     output_lig_mol=%r. '
                                 'Skipping docking to this candidate.' %(tech_prepped_lig_file_list,
                                                                        tech_prepped_prot_file_list,
                                                                        output_receptor_pdb,
                                                                        output_lig_mol))
                    continue
                # Ensure that correct output files exist
                if not(os.path.exists(output_receptor_pdb)) or (os.path.getsize(output_receptor_pdb)==0):
                    logging.info('Docking did not create receptor pdb file %s given inputs:   '
                                 'tech_prepped_lig_file_list=%r   tech_prepped_prot_file_list=%r    '
                                 'output_receptor_pdb=%r     output_lig_mol=%r. Skipping docking '
                                 'to this candidate.' %(output_receptor_pdb,
                                                        tech_prepped_lig_file_list,
                                                        tech_prepped_prot_file_list,
                                                        output_receptor_pdb,
                                                        output_lig_mol))
                    continue
                if not(os.path.exists(output_lig_mol)) or (os.path.getsize(output_lig_mol)==0):
                    logging.info('Docking did not create ligand mol file %s given inputs:   '
                                 'tech_prepped_lig_file_list=%r   tech_prepped_prot_file_list=%r    '
                                 'output_receptor_pdb=%r     output_lig_mol=%r. Skipping docking '
                                 'to this candidate.' %(output_lig_mol,
                                                        tech_prepped_lig_file_list,
                                                        tech_prepped_prot_file_list,
                                                        output_receptor_pdb,
                                                        output_lig_mol))
                    continue
            
                logging.info('Docking was successful for %s. Final receptor and ligand '
                             'files %s and %s exist and are nonzero size.' %(cand_file,
                                                                             output_receptor_pdb,
                                                                             output_lig_mol))

                # Prepare to copy docking results into final result directory
                abs_output_receptor_pdb = os.path.abspath(output_receptor_pdb)
                abs_output_lig_mol = os.path.abspath(output_lig_mol)
            
                # Copy the files one directory up
                os.chdir(abs_targ_dock_dir)
                shutil.copyfile(abs_output_receptor_pdb, output_receptor_pdb)
                shutil.copyfile(abs_output_lig_mol, output_lig_mol)
