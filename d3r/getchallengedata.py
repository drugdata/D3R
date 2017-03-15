#!/usr/bin/env python

import logging
import shutil
import os
import time
import tarfile
from d3r.celpp.filetransfer import WebDavFileTransfer


__author__ = 'j5wagner'



def download_tarball(unpack_dir, config, sleep, max_retry = 50):

    os.chdir(unpack_dir)

    f_f_t_obj = WebDavFileTransfer(config)
    f_f_t_obj.connect()
    retries = 0
    while 1:
        #f_f_t_obj.download_file('/dav/celppweekly/challengedata/latest.txt','./latest.txt')
        #f_f_t_obj.download_file(os.path.join(f_f_t_obj.get_remote_challenge_dir(),'latest.txt'),
        #                        './latest.txt')
        f_f_t_obj.download_file('/dav/challengedata/latest.txt','./latest.txt')
        if os.path.getsize('latest.txt') != 0:
            break
        logging.info('Retry # %i: latest.txt not yet available from challengedata directory' %(retries))
        time.sleep(sleep)
        retries += 1
        if retries == max_retry:
            logging.info('Reached max number of retries (%i). Exiting.' %(max_retry))
            return False
        
    chal_tar_name = open('latest.txt').read().strip()
    #box_chal_tar_name = '/dav/celppweekly/challengedata/' + chal_tar_name
    #box_chal_tar_name = os.path.join(f_f_t_obj.get_remote_challenge_dir(),
    #                                 chal_tar_name)
    box_chal_tar_name = '/dav/challengedata/' + chal_tar_name
    logging.info('Downloaded latest.txt. Challenge package name is %s' %(box_chal_tar_name))
    
    retries = 0
    while 1:
        #f_f_t_obj.download_file(box_chal_tar_name,chal_tar_name)
        f_f_t_obj.download_file(box_chal_tar_name,chal_tar_name)
        if os.path.getsize('latest.txt') != 0:
            break
        logging.info('Retry # %i: %s not yet available from challengedata directory' %(retries, box_chal_tar_name))
        time.sleep(sleep)
        retries += 1
        if retries == max_retry:
            logging.info('Reached max number of retries (%i). Exiting.' %(max_retry))
            return False
        
    return chal_tar_name


def main_get_challenge_data(unpack_dir, ftp_config, local_data_file, sleep):
    abs_orig_dir = os.getcwd()
    abs_unpack_dir = os.path.abspath(unpack_dir)
    if ftp_config is None:
        abs_ftp_config = None
    else:
        abs_ftp_config = os.path.abspath(ftp_config)
    ## Download most recent tarball
    
    if local_data_file is None:
        chal_tar_name = download_tarball(abs_unpack_dir, abs_ftp_config, sleep)
        if chal_tar_name == False:
            logging.info('Unable to download challenge package. Exiting.')
            return False
    else:
        abs_local_data_file = os.path.abspath(local_data_file)
        if not(os.path.exists(abs_local_data_file)):
            logging.info('Specified local data package %s does not exist. Exiting.' %(local_data_file))
            return False
        chal_tar_name = abs_local_data_file

    ## Unpack it in the specified directory
    os.chdir(abs_unpack_dir)
    tarfile_obj = tarfile.open(chal_tar_name, 'r:gz')
    tarfile_obj.extractall()
    tarfile_obj.close()
    return True


if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-u", "--unpackdir", metavar="PATH", help = "Dir where the newest challenge data package should be unpacked")
    parser.add_argument("-f", "--ftpconfig", metavar="PATH", help = "File containing user ftp config information (see included example ftp config for specifics)")
    parser.add_argument("-l", "--localdata", metavar="PATH", help = "Use a local challengedata tarball instead of downloading")
    parser.add_argument("-s", "--sleep", metavar="INT", help = "Number of seconds to wait between checking if challenge package exists", default=3600)

    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level = logging.INFO )
    args = parser.parse_args()
    unpack_dir = args.unpackdir
    ftp_config = args.ftpconfig
    local_data = args.localdata
    sleep = args.sleep
    
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(unpack_dir), 'final.log')

    main_get_challenge_data(unpack_dir, ftp_config, local_data, sleep)

    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)
