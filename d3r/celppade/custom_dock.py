#!/usr/bin/env python

__author__ = 'j5wagner'

import os
import shutil
import re
import glob
import logging

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

    def lig_technical_prep(self, sci_prepped_lig):
        """Technical preparation" is the step immediate preceding
        docking. During this step, you should perform any file
        conversions or processing that are specific to your docking
        program.
        :param sci_prepped_lig: Scientifically prepared ligand file
        :returns: This implementation merely returns the value of
        `sci_prepped_lig` in a list
        """
        return [sci_prepped_lig]

    def receptor_technical_prep(self, sci_prepped_receptor, pocket_center):
        """Technical preparation" is the step immediately preceding
        docking. During this step, you should perform any file
        conversions or processing that are specific to your docking
        program.
        :param sci_prepped_receptor: String containing file name of scientifically prepared receptor from previous step, according to requested suffix.
        :param pocket_center: String containing predicted pocket center.
        :returns: This implementation merely returns the value of 'sci_prepped_receptor" in a list
        """

        # Finally, we return the filenames that will be needed in the
        # docking step. This list is passed to the dock() function as the
        # tech_prepped_receptor_list argument. Here we pass the docking
        # box file (for the docking) and the original scientifically
        # prepped ligand pdb, as that's the easiest way to return the
        # final receptor conformation.
        return [sci_prepped_receptor]

    def dock(self, tech_prepped_lig_list, tech_prepped_receptor_list, output_receptor_pdb, output_lig_mol):
        """# The dock step needs to run the actual docking algorithm. Its first two
        # arguments are the return values from the technical preparation
        # functions for the ligand and receptor. The outputs from this
        # step must be two files - a pdb with the filename specified in
        # the output_receptor_pdb argument, and a mol with the filename
        # specified in the output_ligand_mol argument.
        :returns: Always returns False
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
        pocket_center = [float(i.strip()) for i in pocket_center]
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
        sci_prepped_lig_basename = os.path.basename(sci_prepped_lig_file)
        lig_re_pattern = 'lig_([a-zA-Z0-9]{3})%s' %(Dock.SCI_PREPPED_LIG_SUFFIX)
        lig_re_results = re.findall(lig_re_pattern, sci_prepped_lig_basename)
        if len(lig_re_results) != 1:
            logging.info('Unable to parse prepared ligand %s. Regular expression matching yielded %r' %(sci_prepped_lig_file, lig_re_results))
            return False
        lig_prefix = lig_re_results[0]
        return lig_prefix

    def parse_cand_name(self, cand_name):
        cand_re_pattern = '([a-zA-Z0-9]+)-([a-zA-Z0-9]{4})_([a-zA-Z0-9]{4})%s' %(Dock.SCI_PREPPED_PROT_SUFFIX)
        cand_re_results = re.findall(cand_re_pattern, cand_name)
        if len(cand_re_results) != 1:
            logging.info('Unable to parse prepared protein %s. Regular expression matching yielded %r' %(cand_name, cand_re_results))
            raise Exception('Unable to parse prepared protein %s. Regular expression matching yielded %r' %(cand_name, cand_re_results))
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

            # Get the cand protein names in this directory
            potential_cand_proteins = glob.glob('%s/*-????_????%s' %(targ_prot_prep_dir, Dock.SCI_PREPPED_PROT_SUFFIX))

            # Process potential cand protein names
            for potential_cand_protein in potential_cand_proteins:
                potential_cand_basename = os.path.basename(potential_cand_protein)
                category, targ_id, cand_id = self.parse_cand_name(potential_cand_basename)
            
                # If this is the first valid cand for this target, make a new dictionary key and mark the target as valid
                if not('valid_cands' in targ_dic[targ_name].keys()):
                    targ_dic[targ_name]['valid_targ'] = True
                    targ_dic[targ_name]['prot_prep_dir'] = targ_prot_prep_dir
                    targ_dic[targ_name]['lig_prep_dir'] = targ_lig_prep_dir
                    targ_dic[targ_name]['pocket_center'] = pocket_center
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
            tech_prepped_lig_file_list = self.lig_technical_prep(lig_base_filename)

        
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
                tech_prepped_prot_file_list = self.receptor_technical_prep(prot_base_filename, pocket_center)

                ## Ensure that receptor technical prep was successful
                # Check for function-reported failure
                if tech_prepped_prot_file_list == False:
                    logging.info('Technical protein preparation failed on %s. Skipping candidate %s for target %s.' %(os.path.abspath(prot_base_filename), cand_id, targ_name))
                    continue

                # Ensure that receptor technical prep returns a list of filenames
                if not(type(tech_prepped_prot_file_list) is list):
                    logging.info('Technical protein preparation for %s did not return a list of filenames. Skipping candidate %s for target %s.' %(os.path.abspath(lig_base_filename), cand_id, targ_name))
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
                dock_results = self.dock(tech_prepped_lig_file_list,
                                    tech_prepped_prot_file_list,
                                    output_receptor_pdb,
                                    output_lig_mol)
        
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
