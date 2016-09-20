#!/usr/bin/env python

__author__ = 'j5wagner'

import commands
import os
import glob
import logging
import time
import re 

from d3r.celppade.customeprep import Prep

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.DEBUG )


class OmegaPrep(Prep):

    # This prep script will be required to output files with the appropriate suffixes
    OUTPUT_PROTEIN_SUFFIX = '.pdb'
    OUTPUT_LIG_SUFFIX = '.sdf'



    def ligand_prepare(self, lig_smi_file, out_lig_file, info_dic={}):

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

    OmegaPrep.proteinprep(challenge_data_path, pdb_location, prep_result_path)

    #move the final log file to the result dir
    commands.getoutput("mv %s %s"%(log_file_path, prep_result_path))

