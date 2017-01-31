#!/usr/bin/env python

import logging
import shutil
import os
import glob
import tarfile
from d3r.celpp.filetransfer import WebDavFileTransfer
from d3r.utilities.challenge_data import ChallengeData

__author__ = 'j5wagner'
    

def find_uploadable_results(target_dir):
    # Search through a dock_dir/target_name dir to find files that
    # fit the upload format requirements.
    
    valid_results = []
    abs_targ_dir = os.path.abspath(target_dir)
    potential_pdbs = glob.glob('%s/*-????_????_docked.pdb'%(abs_targ_dir))
    
    for potential_pdb in potential_pdbs:

        ## Do validity checks
        if os.path.getsize(potential_pdb) == 0:
            logging.info('Size of %s is 0. Skipping candidate.'%(potential_pdb))
            continue
        potential_mol = potential_pdb.replace('_docked.pdb','_docked.mol')
        if not(os.path.exists(potential_mol)):
            logging.info('For candidate %s, I was unable to find docked ligand %s Skipping candidate.' %(potential_pdb, potential_mol))
            continue
        if os.path.getsize(potential_mol) == 0:
            logging.info('Size of %s is 0. Skipping candidate.'%(potential_mol))
            continue
        
        # If the target is valid, add it to the packing list
        valid_results.append((potential_pdb, potential_mol))
    
    return valid_results

def make_result_dictionary(dock_dir):
    ## Check every dir that's names like a 4-character target id
    pot_targ_dirs = glob.glob(os.path.join(dock_dir,'????/'))
    result_dic = {}
    for pot_targ_dir in pot_targ_dirs:
        pot_targ_id = os.path.basename(pot_targ_dir.strip('/'))
        valid_results = find_uploadable_results(pot_targ_dir)
        if valid_results == []:
            logging.info('No valid results found for target %s. Skipping target.' %(pot_targ_id))
            continue
        logging.info('Found valid results %r for target %s.' %(valid_results, pot_targ_id))
        result_dic[pot_targ_id] = valid_results
    return result_dic


def main_pack_dock_results(challenge_dir, dock_dir, pack_dir, ftp_config):
    abs_orig_dir = os.getcwd()
    abs_challenge_dir = os.path.abspath(challenge_dir)
    abs_pack_dir = os.path.abspath(pack_dir)

    ## Parse FTP config if given
    if ftp_config is None:
        abs_ftp_config = None
        contestant_id = 'XXXXX'
    else:
        abs_ftp_config = os.path.abspath(ftp_config)
        from d3r.celpp import filetransfer
        f_f_t_obj = WebDavFileTransfer(abs_ftp_config)
        contestant_id = f_f_t_obj.get_contestant_id()

    ## Get this week's name
    chal_data_obj = ChallengeData(abs_challenge_dir)
    week_name = chal_data_obj.get_week_names()[0]

    ## Find all possible uploadable target dirs
    result_dic = make_result_dictionary(dock_dir)
    
    submission_base_name = week_name + '_dockedresults_' + contestant_id 

    
    ## Copy them into this directory
    os.chdir(abs_pack_dir)
    os.mkdir(submission_base_name)
    os.chdir(submission_base_name)
    for targ_id in result_dic:
        os.mkdir(targ_id)
        for docked_pdb, docked_mol in result_dic[targ_id]:
            d_f_basename = os.path.basename(docked_pdb)
            destination = os.path.join(abs_pack_dir,
                                       submission_base_name,
                                       targ_id,
                                       d_f_basename)
            shutil.copyfile(docked_pdb, destination)
            
            d_f_basename = os.path.basename(docked_mol)
            destination = os.path.join(abs_pack_dir,
                                       submission_base_name,
                                       targ_id,
                                       d_f_basename)
            shutil.copyfile(docked_mol, destination)
            
            ### TEMPORARY HACK UNTIL EVALUATE GETS FIXED ###
            #hack_folder_name = destination.replace('_docked.mol','')
            #os.mkdir(hack_folder_name)
            #destination = os.path.join(hack_folder_name, d_f_basename)
            #shutil.copyfile(docked_mol, destination)
            #destination = destination.replace('.mol','.pdb')
            #shutil.copyfile(docked_pdb, destination)

    
    
    ## Tar up the pack directory. To keep this tarball from containing
    ## absolute paths, we enter the pack directory, go one directory
    ## up, and pack the tarball from there (thereby giving the content
    ## reasonable relative paths)
    os.chdir(abs_pack_dir)
    #os.chdir('..')
    
    #abs_tar_name = abs_pack_dir.rstrip('/') + '.tar.gz'
    tar_base_name = submission_base_name + '.tar.gz'
    abs_tar_name = os.path.join(abs_pack_dir, tar_base_name)
    
    logging.info('Creating tarfile %s in directory %s' %(abs_tar_name, os.getcwd()))
    tarfile_obj = tarfile.open(abs_tar_name, 'w:gz')
    logging.info('Writing to tarfile')
    #tarfile_obj.add(os.path.basename(abs_pack_dir.rstrip('/')))
    tarfile_obj.add(submission_base_name)
    tarfile_obj.close()
    logging.info('Tarfile closed')


    ## Use ftp config to upload tarball
    if ftp_config is None:
        logging.info('No ftp_config file given. Skipping upload')
        return

    #tar_base_name = os.path.basename(abs_tar_name)
    f_f_t_obj.connect()
    submission_dir = os.path.join('/dav/',
                                  #f_f_t_obj.get_remove_submission_dir, 
                                  contestant_id)
    f_f_t_obj.upload_file_direct(abs_tar_name,
                                 submission_dir,
                                 #f_f_t_obj.get_remote_submission_dir(),
                                 #'/dav/celppweekly/usersubmissions/12345/',
                                 tar_base_name)
    logging.info(f_f_t_obj.get_upload_summary())
    f_f_t_obj.disconnect()


        

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-c", "--challengedata", metavar="PATH", help = "PATH to the unpacked challenge data package")
    parser.add_argument("-d", "--dockdir", metavar = "PATH", help = "Dir where docking was performed")
    parser.add_argument("-p", "--packdir", metavar="PATH", help = "Dir where the packing and uploading will be performed. Note that, for competition entries, this expects this dir name to already have the proper formatting for this contestant's entry, and will be used as the base of the tar.gz file.")
    parser.add_argument("-f", "--ftpconfig", metavar="PATH", help = "File containing user ftp config information (see included example ftp config for specifics)")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level = logging.INFO )
    args = parser.parse_args()
    challenge_dir = args.challengedata
    dock_dir = args.dockdir
    pack_dir = args.packdir
    ftp_config = args.ftpconfig
    
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(pack_dir), 'final.log')

    main_pack_dock_results(challenge_dir, dock_dir, pack_dir, ftp_config)

    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)
