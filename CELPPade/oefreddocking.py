#!/usr/bin/env python

import commands
import os
import shutil
import re
import glob
import logging
import time
from d3r.celpp import task

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )

################################################################
########## CONTESTANT MODIFICATIONS SHOULD BEGIN HERE ##########
################################################################

## To ensure that the workflow copies over the correct scientific prep
## outputs, list the expected filetypes here. This script will look
## for files in the target prep directory that follow the naming
## convention "<cand category>-<target pdbid>_<cand
## pdbid>_prepared.<sci_prepped_prot_suffix>", for example
## "LMCSS-5b4t_3w8e_prepared.pdb". It will then pass the appropriate
## filenames and prefixes to the technical_prep() and dock()
## functions.
sci_prepped_lig_suffix = '_prepared.sdf'
sci_prepped_prot_suffix = '_prepared.pdb'

def lig_technical_prep(sci_prepped_lig):
    # "Technical preparation" is the step immediate preceding
    # docking. During this step, you should perform any file
    # conversions or processing that are specific to your docking
    # program.

    # For FRED, there is no processing needed for the ligand beyond
    # conformer generation. In our workflow, this was performed in the
    # scientific prep step. 

    # Return values from this function must always be lists, even if
    # there is only one file.
    return [sci_prepped_lig]

def receptor_technical_prep(sci_prepped_receptor, pocket_center):
    # "Technical preparation" is the step immediately preceding
    # docking. During this step, you should perform any file
    # conversions or processing that are specific to your docking
    # program.

    # For FRED, receptor_technical_prep takes a scientifically
    # prepared pdb file and creates a docking grid.

    receptor_prefix = sci_prepped_receptor.replace('.pdb','')
    docking_box_file = receptor_prefix + '.oeb.gz'

    # It helps to capture stdout and stderr from these steps for debugging
    stdout_file = receptor_prefix + '_docking_box_gen_stdout'
    stderr_file = receptor_prefix + '_docking_box_gen_stderr'

    # Technical prep for FRED consists of generating a docking box. The
    # box is defined by a .xyz file with two atoms specifying the
    # corners
    box_xyz_file = 'box.xyz'
    center_list = pocket_center.split(',')
    box_size = 30.
    corner_1 = (float(center_list[0]) - (box_size / 2),
                float(center_list[1]) - (box_size / 2),
                float(center_list[2]) - (box_size / 2))
    corner_2 = (float(center_list[0]) + (box_size / 2),
                float(center_list[1]) + (box_size / 2),
                float(center_list[2]) + (box_size / 2))
    box_xyz_text = '''2 
CELPP auto-generated docking box 
 X %i %i %i 
 X %i %i %i 
''' %(corner_1[0], corner_1[1], corner_1[2], 
      corner_2[0], corner_2[1], corner_2[2])
    with open(box_xyz_file, 'wb') as of:
        of.write(box_xyz_text)
        
    # Run the box generation command
    cmd = 'receptor_setup -protein ' + sci_prepped_receptor + ' -box ' + box_xyz_file + ' -receptor ' + docking_box_file + ' 2> ' + stderr_file + ' 1> ' + stdout_file
    commands.getoutput(cmd)
    
    # Check for succeess
    if not(os.path.exists(docking_box_file)):
        # This step will be run on all ligands, so it helps to be very descriptive about file locations
        abs_stdout_file = os.path.abspath(stdout_file)
        abs_stderr_file = os.path.abspath(stderr_file)
        logging.info('Docking box generation failed. See %s and %s for details' %(abs_stdout_file, abs_stderr_file))
        # CELPPade handles errors by returning False
        return False

    # Finally, we return the filenames that will be needed in the
    # docking step. This list is passed to the dock() function as the
    # tech_prepped_receptor_list argument. Here we pass the docking
    # box file (for the docking) and the original scientifically
    # prepped ligand pdb, as that's the easiest way to return the
    # final receptor conformation.
    return [docking_box_file, sci_prepped_receptor]

    


def dock(tech_prepped_lig_list, tech_prepped_receptor_list, output_receptor_pdb, output_lig_mol):
    # The dock step runs the actual docking algorithm. Its first two
    # arguments are the return values from the technical preparation
    # functions for the ligand and receptor. The outputs from this
    # step must be two files - a pdb with the filename specified in
    # the output_receptor_pdb argument, and a mol with the filename
    # specified in the output_ligand_mol argument.

    tech_prepped_grid = tech_prepped_receptor_list[0]
    sci_prepped_receptor = tech_prepped_receptor_list[1]
    receptor_prefix = os.path.basename(tech_prepped_grid).replace('.oeb.gz','')

    tech_prepped_lig = tech_prepped_lig_list[0]

    # FRED can only output .oeb and .sdf files. We'll have it output a
    # .sdf, then convert it to mol after FRED runs.
    docked_lig_sdf = receptor_prefix + '_docked.sdf'
    #docked_lig_mol = receptor_prefix + '_docked.mol'

    # It helps to capture stdout and stderr from these steps for debugging
    fred_stdout_file = receptor_prefix + '_fred_dock_stdout'
    fred_stderr_file = receptor_prefix + '_fred_dock_stderr'
    babel_stdout_file = receptor_prefix + '_babel_sdf2mol_stdout'
    babel_stderr_file = receptor_prefix + '_babel_sdf2mol_stderr'


    # Run the docking command
    dock_command = 'fred -receptor ' + tech_prepped_grid + ' -dbase ' + tech_prepped_lig + ' -docked_molecule_file ' + docked_lig_sdf + ' 2> ' + fred_stderr_file + ' 1> ' + fred_stdout_file
    commands.getoutput(dock_command)
    if not(os.path.isfile(docked_lig_sdf)):
        logging.info("Docking failed for protein %s and ligand %s" %(tech_prepped_grid, tech_prepped_lig))
        return False
    
    ## Create the required output files
    # First, the ligand needs to be converted from .sdf to .mol
    lig_convert_command = 'babel -isdf ' + docked_lig_sdf + ' -omol ' + output_lig_mol +  ' 2> ' + babel_stderr_file + ' 1> ' + babel_stdout_file
    commands.getoutput(lig_convert_command)
    if not(os.path.isfile(output_lig_mol)):
        logging.info("Babel conversion failed for docked ligand  %s" %(docked_lig_sdf))
        return False

    # And now we need to prepare the final receptor file. Because this
    # docking didn't modify the receptor, we will just copy the
    # original scientifically-prepped receptor pdb and give it the
    # appropriate name.
    shutil.copyfile(sci_prepped_receptor, output_receptor_pdb)

    # And do a final sanity check
    if not(os.path.isfile(output_receptor_pdb)):
        logging.info("File copy failed for receptor %s" %(output_receptor_pdb))
        return False
    
    return True

################################################################
########### CONTESTANT MODIFICATIONS SHOULD END HERE ###########
################################################################

def get_pocket_center(target_prep_dir):
    pocket_center_file = os.path.join(target_prep_dir,'center.txt')
    try:
        pocket_center = open(pocket_center_file,"r").readlines()[0]
    except:
        logging.info("Unable to find the center file for this case %s."%targ_name)
        return False
    return pocket_center

def get_sci_prepped_lig(target_prep_dir, sci_prepped_lig_suffix):
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


def parse_lig_filename(sci_prepped_lig_file):
    sci_prepped_lig_basename = os.path.basename(sci_prepped_lig_file)
    lig_re_pattern = 'lig_([a-zA-Z0-9]{3})%s' %(sci_prepped_lig_suffix)
    lig_re_results = re.findall(lig_re_pattern, sci_prepped_lig_basename)
    if len(lig_re_results) != 1:
        logging.info('Unable to parse prepared ligand %s. Regular expression matching yielded %r' %(lig_re_results))
        return False
    lig_prefix = lig_re_results[0]
    return lig_prefix


def parse_cand_name(cand_name):
    cand_re_pattern = '([a-zA-Z0-9]+)-([a-zA-Z0-9]{4})_([a-zA-Z0-9]{4})%s' %(sci_prepped_prot_suffix)
    cand_re_results = re.findall(cand_re_pattern, cand_name)
    if len(cand_re_results) != 1:
        logging.info('Unable to parse prepared protein %s. Regular expression matching yielded %r' %(cand_name, cand_re_results))
        raise Exception('Unable to parse prepared protein %s. Regular expression matching yielded %r' %(cand_name, cand_re_results))
    category, target, cand = cand_re_results[0]
    return category, target, cand



def main_dock (prep_result_dir, dock_dir, update= True):
    #os.chdir(prep_result_dir)
    abs_prep_dir = os.path.abspath(prep_result_dir)
    abs_dock_dir = os.path.abspath(dock_dir)

    targ_prep_dirs = glob.glob('%s/????' %(abs_prep_dir))
    targ_dic = {}
    
    for targ_prep_dir in targ_prep_dirs:
        targ_name = os.path.basename(targ_prep_dir.strip('/'))
        targ_dock_dir = os.path.join(abs_dock_dir,targ_name)
        #os.mkdir(target_dock_dir)
        logging.info("============= Starting to process target:%s =============" %targ_name)
        targ_dic[targ_name] = {}
        # Until we find a good ligand and cand receptor, mark this as invalid
        targ_dic[targ_name]['valid_targ'] = False
        
        # Get the binding pocket center
        pocket_center = get_pocket_center(targ_prep_dir)
        if pocket_center == False:
            logging.info('Failed to find pocket center file in dirctory %s. Skipping target %s.' %(targ_prep_dir, targ_name))
            continue

        # Get the ligand name in this directory
        sci_prepped_lig_file = get_sci_prepped_lig(targ_prep_dir, 
                                                      sci_prepped_lig_suffix)
        if sci_prepped_lig_file == False:
            logging.info('Unable to find single ligand for target dir %s. Skipping target %s.' %(targ_prep_dir, targ_name))
            continue

        # Process ligand file names
        lig_prefix = parse_lig_filename(sci_prepped_lig_file)
        if lig_prefix == False:
            logging.info('Unable to parse ligand filename %s. Skipping target %s.' %(sci_prepped_lig_file, targ_name))

        # Get the cand protein names in this directory
        potential_cand_proteins = glob.glob('%s/*-????_????%s' %(targ_prep_dir, sci_prepped_prot_suffix))

        # Process potential cand protein names 
        for potential_cand_protein in potential_cand_proteins:
            potential_cand_basename = os.path.basename(potential_cand_protein)
            category, targ_id, cand_id = parse_cand_name(potential_cand_basename)
            
            # If this is the first valid cand for this target, make a new dictionary key and mark the target as valid
            if not('valid_cands' in targ_dic[targ_name].keys()):
                targ_dic[targ_name]['valid_targ'] = True
                targ_dic[targ_name]['prep_dir'] = targ_prep_dir
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
        tech_prepped_lig_file_list = lig_technical_prep(lig_base_filename)

        
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
            tech_prepped_prot_file_list = receptor_technical_prep(prot_base_filename, pocket_center)

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
            dock_results = dock(tech_prepped_lig_file_list,
                                tech_prepped_prot_file_list,
                                output_receptor_pdb,
                                output_lig_mol)
        
            ## Check for success
            # Check for self-reported failure
            if dock_results == False:
                logging.info('Docking returned False given inputs: tech_prepped_lig_file_list=%r   tech_prepped_prot_file_list=%r    output_receptor_pdb=%r     output_lig_mol=%r. Skipping docking to this candidate.' %(tech_prepped_lig_file_list,
                                                                                                                                                                                  tech_prepped_prot_file_list,
                                                                                                                                                                                  output_receptor_pdb,
                                                                                                                                                                                  output_lig_mol))
                continue
            # Ensure that correct output files exist
            if not(os.path.exists(output_receptor_pdb)):
                logging.info('Docking did not create receptor pdb file %s given inputs:'
                             'tech_prepped_lig_file_list=%r   tech_prepped_prot_file_list=%r    '
                             'output_receptor_pdb=%r     output_lig_mol=%r. Skipping docking '
                             'to this candidate.' %(output_receptor_pdb,
                                                    tech_prepped_lig_file_list,
                                                    tech_prepped_prot_file_list,
                                                    output_receptor_pdb,
                                                    output_lig_mol))
                continue
            if not(os.path.exists(output_lig_mol)):
                logging.info('Docking did not create ligand mol file %s given inputs:'
                             'tech_prepped_lig_file_list=%r   tech_prepped_prot_file_list=%r    '
                             'output_receptor_pdb=%r     output_lig_mol=%r. Skipping docking '
                             'to this candidate.' %(output_lig_mol, 
                                                    tech_prepped_lig_file_list,
                                                    tech_prepped_prot_file_list,
                                                    output_receptor_pdb,
                                                    output_lig_mol))
                continue
            
            logging.info('Docking was successful for %s. Final receptor and ligand files %s and %s exist' %(cand_file, output_receptor_pdb, output_lig_mol))
            # Prepare to copy docking results into final result directory
            abs_output_receptor_pdb = os.path.abspath(output_receptor_pdb)
            abs_output_lig_mol = os.path.abspath(output_lig_mol)
            
            # Copy the files one directory up
            os.chdir(abs_targ_dock_dir)
            shutil.copyfile(abs_output_receptor_pdb, output_receptor_pdb)
            shutil.copyfile(abs_output_lig_mol, output_lig_mol)


'''
        for cand_protein in cand_proteins:
            logging.info( "Working on this receptor: %s"%cand_protein )
            if not os.path.isdir(cand_prefix):
                os.mkdir(cand_prefix)
            os.chdir(cand_prefix)
            commands.getoutput("cp ../%s ."%cand_protein)
            
            ## Technical prep: Prepare the protein
            #receptorMol2File = '.'.join(possible_protein.split('.')[:-1]) + '.mol2'
            commands.getoutput('. /usr/local/mgltools/bin/mglenv.sh; pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r %s' %(cand_protein))
            receptorPdbqtFile = cand_protein.replace('.mol2','.pdbqt')


            for ligand_mol2 in ligand_mol2s:
                ## Technical prep: Prepare the ligand
                ligand_name = ligand_mol2.replace('lig_','').replace('_prepped.mol2','')
                ligandPdbqtFile = ligand_mol2.replace('.mol2','.pdbqt')
                commands.getoutput("cp ../%s ." %(ligand_mol2))
                commands.getoutput('. /usr/local/mgltools/bin/mglenv.sh; pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l %s' %(ligand_mol2))
    
                ## Do the docking
                logging.info("Trying to dock...")
                out_dock_file = dock(ligandPdbqtFile, receptorPdbqtFile, grid_center)
                if out_dock_file == False:
                    logging.info("Docking failed - Skipping")
                    continue
                logging.info("Finished docking, beginning post-docking file conversion")
                
                ## Convert the receptor to pdb
                intermediates_prefix = '%s_postdocking' %(cand_prefix)
                output_prefix = '%s_docked' %(cand_prefix)
                #receptorPdbqt = out_dock_file.replace('ligand_out',cand_name)
                #receptorPdbqt = cand_prefix+'.pdbqt'
                ## This receptor pdb will be one of our final outputs
                outputReceptorPdb = "%s.pdb" %(output_prefix)
                commands.getoutput('. /usr/local/mgltools/bin/mglenv.sh; python $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f %s -o %s' %(receptorPdbqtFile, outputReceptorPdb))
                
                ## Then make the ligand mol
                ## pdbqt_to_pdb.py can't split up the multiple poses in vina's output files, so we do that by hand
                with open(out_dock_file) as fo:
                    fileData = fo.read()
                fileDataSp = fileData.split('ENDMDL')
                
                ## Write out each pose to its own pdb and mol file, then merge with the receptor to make the complex files.
                for index, poseText in enumerate(fileDataSp[:-1]):
                    this_pose_pdbqt = intermediates_prefix+'_ligand'+str(index+1)+'.pdbqt'
                    this_pose_pdb = intermediates_prefix+'_ligand'+str(index+1)+'.pdb'
                    this_pose_mol = intermediates_prefix+'_ligand'+str(index+1)+'.mol'
                    with open(this_pose_pdbqt,'wb') as wf:
                        wf.write(poseText+'ENDMDL')
                    commands.getoutput('. /usr/local/mgltools/bin/mglenv.sh; python $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f %s -o %s' %(this_pose_pdbqt, this_pose_pdb))
                    commands.getoutput('babel -ipdb %s -omol %s' %(this_pose_pdb, this_pose_mol))
                
                ## Here convert our top-ranked pose to the final submission for this docking
                ## Right now we're ignoring everything other than the top pose
                top_intermediate_mol = intermediates_prefix+'_ligand1.mol'
                final_ligand_mol = output_prefix+'.mol'
                commands.getoutput('cp %s %s' %( top_intermediate_mol, final_ligand_mol))
                    
                ##################
            os.chdir(abs_target_dir)
        os.chdir(abs_dock_dir)
        ##################
'''

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-s", "--sciprepdir", metavar="PATH", help = "PATH where we can find the scientific protein and lig prep output")
    parser.add_argument("-o", "--outdir", metavar = "PATH", help = "PATH where we will put the docking output")
    parser.add_argument("-u", "--update", action = "store_true",  help = "Update the docking result", default = False) 
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )
    opt = parser.parse_args()
    sci_prep_dir = opt.sciprepdir
    dock_dir = opt.outdir
    update = opt.update
    #running under this dir
    running_dir = os.getcwd()
    main_dock(sci_prep_dir, dock_dir, update = update)
    #move the final log file to the result dir
    #log_file_path = os.path.join(running_dir, 'final.log')
    #commands.getoutput("mv %s %s"%(log_file_path, dock_dir))
