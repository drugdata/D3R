#!/usr/bin/env python

import commands
import os
import glob
import logging
import time

logger = logging.getLogger()
logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )

pdb2MaePyText = '''import schrodinger.structure as ss
import sys
inFile = sys.argv[1]
outFile = sys.argv[2]
if not('.mae' in inFile):
    raise Exception('Invalid input file; must be .mae format')
struct = ss.StructureReader(inFile).next()
print struct
struct.write(outFile)
'''

def dock(ligand_pdbqt,  protein_pdbqt, grid_center,):
    center_list = grid_center.split(',')
    vina_command = 'vina --receptor %s  --ligand %s --center_x %s --center_y %s --center_z %s --size_x 15 --size_y 15 --size_z 15 --seed 999 >& vina_output' %(protein_pdbqt, ligand_pdbqt, center_list[0], center_list[1], center_list[2])
    out_dock_file = ligand_pdbqt.replace('.pdbqt','_out.pdbqt')

    commands.getoutput(vina_command)
    if not(os.path.isfile(out_dock_file)):
        logging.info("Docking failed for protein %s and ligand %s" %(protein_pdbqt, ligand_pdbqt))
        return False
    return out_dock_file

def main_vina (stage_3_result, stage_4_working, update= True):
    os.chdir(stage_4_working)
    current_dir = os.getcwd()
    dockable_paths = []
    for all_pdb_path in os.walk(stage_3_result):
        if all_pdb_path[0] != stage_3_result:
            dockable_paths.append(all_pdb_path[0])

    for dockable_path in dockable_paths:
        commands.getoutput("cp -r %s %s" %(dockable_path, stage_4_working))
        target_name = os.path.basename(dockable_path)
        os.chdir(target_name)
        current_dir_layer_2 = os.getcwd()
        ##################
        #Do the docking here
        
        #1, get the center
        try:
            grid_center = open("center.txt","r").readlines()[0]
            logging.info("=============Start working on this case:%s=========="%target_name)
        except:
            logging.info("Fatal error: Unable to find the center file for this case %s"%target_name)
            os.chdir(current_dir)
            continue
        
        # Get the candidate protein names in this directory
        candidate_proteins = glob.glob('./*-????_????_prepared.mol2')
        
        logging.info('Found candidates %r' %(candidate_proteins))
        # Get the ligand names in this directory
        ligand_mol2s = glob.glob('lig_*_prepared.mol2')
        ligand_mol2s = [i for i in ligand_mol2s if not "unprep" in i]
        if len(ligand_mol2s) == 0:
            logging.info('No ligand files found for target %s. Skipping.' %(dockable_path))
            os.chdir(current_dir)
            continue
        if len(ligand_mol2s) > 1:
            logging.info('Multiple ligand files found for target %s (found ligands %r). The workflow currently should only be sending one ligand. Skipping. ' %(dockable_path, ligand_mol2s))
            os.chdir(current_dir)
            continue

        
        
        for candidate_protein in candidate_proteins:
            candidate_prefix = candidate_protein.replace('_prepared.mol2','')
            #if os.path.isfile(candidate_protein):
            logging.info( "Working on this receptor: %s"%candidate_protein )
            if not os.path.isdir(candidate_prefix):
                os.mkdir(candidate_prefix)
            os.chdir(candidate_prefix)
            ##################
            #do docking inside
            commands.getoutput("cp ../%s ."%candidate_protein)
            
            ## Technical prep: Prepare the protein
            #receptorMol2File = '.'.join(possible_protein.split('.')[:-1]) + '.mol2'
            commands.getoutput('. /usr/local/mgltools/bin/mglenv.sh; pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r %s' %(candidate_protein))
            receptorPdbqtFile = candidate_protein.replace('.mol2','.pdbqt')


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
                intermediates_prefix = '%s_postdocking' %(candidate_prefix)
                output_prefix = '%s_docked' %(candidate_prefix)
                #receptorPdbqt = out_dock_file.replace('ligand_out',candidate_name)
                #receptorPdbqt = candidate_prefix+'.pdbqt'
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
                commands.getoutput('cp %s ../' %(final_ligand_mol))
                commands.getoutput('cp %s ../' %(outputReceptorPdb))
                    
                ##################
            os.chdir(current_dir_layer_2)
        os.chdir(current_dir)
        ##################

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-s", "--structuredir", metavar="PATH", help = "PATH where we can find the proteinligprep output")
    parser.add_argument("-o", "--outdir", metavar = "PATH", help = "PATH where we will put the docking output")
    parser.add_argument("-u", "--update", action = "store_true",  help = "Update the docking result", default = False) 
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )
    opt = parser.parse_args()
    proteinligprep_dir = opt.structuredir
    dock_result_dir = opt.outdir
    update = opt.update
    #running under this dir
    running_dir = os.getcwd()
    main_vina(proteinligprep_dir, dock_result_dir, update = update)
    #move the final log file to the result dir
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s"%(log_file_path, dock_result_dir))
