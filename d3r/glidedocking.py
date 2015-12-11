#!/usr/bin/evn python

import commands
import os
import glob
import logging
import time

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )


def grid (grid_center, prep_pro):
    protein_name = prep_pro.split(".")[0]
    out_grid = protein_name + ".zip"
    grid_in  = protein_name +  "_grid.in"
    grid_1 = "GRID_CENTER \t %s \n"%grid_center
    grid_2 = "GRIDFILE \t %s \n"%out_grid
    grid_3 = "INNERBOX \t 10, 10, 10 \n"
    grid_4 = "OUTERBOX \t 30, 30, 30 \n"
    grid_5 = "RECEP_FILE \t %s \n"%prep_pro
    f = open (grid_in, "w")
    for grid_line in (grid_1, grid_2, grid_3, grid_4, grid_5):
        f.writelines(grid_line)
    f.close()
    commands.getoutput("$SCHRODINGER/glide -WAIT %s"%grid_in)
    if os.path.isfile(out_grid):
        return out_grid
    else:
        return False

def dock(ligand_file, out_grid, protein_title):
    dock_in = protein_title + "_dock.in"
    dock_1 = "GRIDFILE \t %s \n"%out_grid
    dock_2 = "LIGANDFILE \t %s \n"%ligand_file
    dock_3 = "POSTDOCK_XP_DELE \t 0.5 \n"
    dock_4 = "PRECISION \t XP \n"
    dock_5 = "EXPANDED_SAMPLING \t True \n"
    dock_6 = "POSES_PER_LIG \t 10 \n"
    dock_7 = "WRITE_XP_DESC \t False \n"
    f = open(dock_in, "w")
    for dock_line in (dock_1, dock_2, dock_3, dock_4, dock_5, dock_6, dock_7):
        f.writelines(dock_line)
    f.close()
    out_dock_file = protein_title + "_dock_pv.maegz"
    commands.getoutput("$SCHRODINGER/glide -WAIT %s"%dock_in)
def main_glide (stage_3_result, stage_4_working):
    os.chdir(stage_4_working)
    current_dir = os.getcwd()
    dockable_paths = []
    for all_pdb_path in os.walk(stage_3_result):
        if all_pdb_path[0] != stage_3_result:
            dockable_paths.append(all_pdb_path[0])

    for dockable_path in dockable_paths:
        commands.getoutput("cp -r %s %s"%(dockable_path, stage_4_working))
        dockable_path_local = os.path.basename(dockable_path)
        os.chdir(dockable_path_local)
        current_dir_layer_2 = os.getcwd()
        ##################
        #Do the docking here
        #1, get the center
        grid_center = open("center.txt","r").readlines()[0]
        if grid_center:
            
            logging.info("=============Start working in this case:%s=========="%dockable_path_local)
        else:
            logging.info("Fatal error: Unable to find the center file for this case %s"%dockable_path_local)
            continue
        for possible_protein in ("largest.maegz", "smallest.maegz", "holo.maegz", "apo.maegz"):
        #for possible_protein in ("smallest.maegz", "holo.maegz", "apo.maegz"):
            if os.path.isfile(possible_protein):
                logging.info( "Working in this pdb: %s"%possible_protein )
                protein_title = possible_protein.split(".")[0]
                if not os.path.isdir(protein_title):
                    os.mkdir(protein_title)
                os.chdir(protein_title)
                ##################
                #do docking inside
                commands.getoutput("cp ../%s ."%possible_protein)
                commands.getoutput("cp ../ligand.mae .")
                #generate grid
                out_grid = grid(grid_center, possible_protein) 
                logging.info("Trying to get the grid file...")
                time.sleep(200)
                logging.info("Finish grid generation...")
                dock("ligand.mae", out_grid, protein_title)
                logging.info("Finish docking...")
                ##################
                os.chdir(current_dir_layer_2)
        os.chdir(current_dir)
        #raw_input()
        ##################

if ("__main__") == (__name__):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--structuredir", metavar="PATH", help = "PATH where we could find the stage 2 output")
    parser.add_option("-o", "--outdir", metavar = "PATH", help = "PATH where we run stage 3")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )
    (opt, args) = parser.parse_args()
    stage_3_result = opt.structuredir
    stage_4_result = opt.outdir
    #running under this dir
    running_dir = os.getcwd()
    main_glide(stage_3_result, stage_4_result)
    #move the final log file to the result dir
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s"%(log_file_path, stage_4_result))
