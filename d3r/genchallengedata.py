#!/usr/bin/env python
__author__ = 'sliu'
#from stage 3 blastnfilter result generating the stage 4 dataset to release to the participant
import os
import sys
import glob
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
import commands


pymol_align_text = '''
import __main__
__main__.pymol_argv = [ 'pymol', '-c'] #put -cq here to suppress pymol output
import pymol
pymol.finish_launching()
import sys

caAlignment = True
proxim_filter_dist = 20

#molPrefixes = ['LMCSS_splitted_receptor1',
#               'SMCSS_splitted_receptor1',
#               'hiResApo_splitted_receptor1',
#               'hiResHolo_splitted_receptor1',]

reference_pdb = sys.argv[1]
mobile_pdb = sys.argv[2]
print 'reference_pdb is', reference_pdb
print 'mobie_pdb is', mobile_pdb

reference_pdb_obj_name = reference_pdb.split('/')[-1].split('.pdb')[0]
mobile_pdb_obj_name = mobile_pdb.split('/')[-1].split('.pdb')[0]

pymol.cmd.load(reference_pdb)
pymol.cmd.load(mobile_pdb)


### Geometric filter ###
# >>> help(pymol.cmd.pseudoatom)
#Help on function pseudoatom in module pymol.creating:
#
#pseudoatom(object='', selection='', name='PS1', resn='PSD', resi='1', chain='P', segi='PSDO', elem='PS', vdw=-1.0, hetatm=1, b=0.0, q=0.0, color='', label='', pos=None, state=0, mode='rms', quiet=1,)

center = open('center.txt').read().split()
center = [float(i.strip(',')) for i in center]
print center
pymol.cmd.pseudoatom('LMCSS_lig_center', pos=center)
pymol.cmd.select('reference',"(/%s////CA) and (byres (LMCSS_lig_center around %i))" %(reference_pdb_obj_name,proxim_filter_dist) )


print
print
print '=========================='
print "PROCESSING", mobile_pdb_obj_name


### Secondary structure filter ###  
# Get resnums and their SS types
pymol.stored.tuples = []
resiSS = pymol.cmd.iterate(mobile_pdb_obj_name,'stored.tuples.append((resn,resv,ss))')

# Filter to "H" and "S" with list comp
print "FILTERING FOR SS"
ss_near_bs = [i for i in pymol.stored.tuples if i[2] in ['H','S']]
resv_sel_str = '+'.join([str(i[1]) for i in ss_near_bs])


### 3D alignment ###


##  API reference  ##
## http://www.pymolwiki.org/index.php/Align
#cmd.align( string mobile, string target, float cutoff=2.0,
#           int cycles=5, float gap=-10.0, float extend=-0.5,
#           int max_gap=50, string object=None, string matrix='BLOSUM62',
#           int mobile_state=0, int target_state=0, int quiet=1,
#           int max_skip=0, int transform=1, int reset=0 )


if caAlignment:
    alnOut = pymol.cmd.align('/%s////CA'%(mobile_pdb_obj_name),
                             'reference'
                             )
else:
    alnOut = pymol.cmd.align('/%s////'%(mobile_pdb_obj_name),
                             'reference',
                             )
pymol.cmd.save('rot-%s.pdb' %(mobile_pdb_obj_name),
               '/%s' %(mobile_pdb_obj_name)
               )
print alnOut
pymol.cmd.save('alignment_reference.pdb' ,
               'reference'
               )

#print sys.stdout.read()
pymol.cmd.quit()
'''

                                                                                                                                                                                           

def parse_txt (txt_filename):
    #parsing the result txt file and give back the

    #1, inchi string for the ligand 
    #2, LMCSS, SMCSS, hiResHolo, hiResApo structure ID
    info_dic = {}
    target_name_txt = os.path.basename(txt_filename).split(".")[0]
    txt_file = open(txt_filename, "r")
    lines = txt_file.readlines()
    txt_file.close()
    for line in lines:
        info_title = line.split(",")[0]
        info_value_1 = line.split(",")[1]
        if line.split(",")[0] == "inchi":
            info_title = line.split(", InChI=")[0]
            info_value_1 = line.split(", InChI=")[1].split("\n")[0]
        else:
            try:
                info_value_2 = line.split(",")[2].split("\n")[0].strip()
                info_title = line.split(",")[0]
                info_value_1 = line.split(",")[1].strip()
            except:
                info_value_2 = False
                info_title = line.split(",")[0]
                info_value_1 = line.split(",")[1].split("\n")[0].strip()
                pass
        if not info_value_2:                                  
            #just use the first item if there are multiple entries (use the highest resolution one)
            if info_title not in info_dic:                    
                info_dic[info_title] = info_value_1           
        else:                                                 
            if info_title not in info_dic:                    
                info_dic[info_title] = (info_value_1,info_value_2)           
    #check if the query protein id is identical with the filename
    if not info_dic["query"] == target_name_txt:
        logging.info("The query pdb ID:%s is different from the filename: %s, need to check"%(info_dic["query"], target_name_txt)) 
    return info_dic  

def align_proteins (target_protein, pre_prepare_protein, post_prepare_protein):
    #commands.getoutput("$SCHRODINGER/utilities/structalign %s %s"%(target_protein, pre_prepare_protein))
    with open('pymolAlign.py','wb') as of:
        of.write(pymol_align_text)
    stdout_file = 'pymol_align_%s_out' %(pre_prepare_protein)
    stderr_file = 'pymol_align_%s_err' %(pre_prepare_protein)
    cmd = "python pymolAlign.py %s %s 2> %s 1> %s"%(target_protein, pre_prepare_protein, stderr_file, stdout_file)
    logging.info('Running alignment command %s' %(cmd))
    commands.getoutput(cmd)

    rotated_protein = "rot-" + pre_prepare_protein
    stderr_output = open(stderr_file).read().strip()
    if os.path.isfile(rotated_protein) and stderr_output== '':
        commands.getoutput("mv %s %s"%(rotated_protein, post_prepare_protein))
        return True
    else:
        logging.info('Alignment failed. Command was %s. Stderr output was %r' %(cmd, stderr_output))
        return False

def generate_ligand (inchi, ligand_title, ):
    valid_inchi = "InChI="+inchi
    rd_mol = Chem.MolFromInchi(format(valid_inchi), removeHs=False, sanitize=False, treatWarningAsError=True)
    smiles = Chem.MolToSmiles(rd_mol, isomericSmiles=True)
    smile_filename = "lig_" + ligand_title + ".smi"
    smile_file = open(smile_filename, "w")
    smile_file.writelines(smiles)
    inchi_filename = "lig_" + ligand_title + ".inchi"
    inchi_file = open(inchi_filename, "w")
    inchi_file.writelines(valid_inchi)
    mol_filename = "lig_" + ligand_title + ".mol"
    #rd_mol_H = Chem.AddHs(rd_mol)
    #AllChem.EmbedMolecule(rd_mol_H)
    #AllChem.UFFOptimizeMolecule(rd_mol_H)
    AllChem.Compute2DCoords(rd_mol)
    Chem.MolToMolFile (rd_mol, mol_filename, includeStereo=True)

def pull_ligand_out (proteinfile, ligname, ligandfile):
    p_xyz_lines = open(proteinfile,"r").readlines()
    multi_ligand = False
    l_xyz_lines = []
    atom_list = []
    for p_xyz_line in p_xyz_lines:
        if "HETATM" in p_xyz_line and ligname in p_xyz_line[17:21]:
            atom_name = p_xyz_line[12:16]
            if atom_name in atom_list:
                multi_ligand = True
            else:
                atom_list.append(atom_name)
                l_xyz_lines.append(p_xyz_line)
    if not multi_ligand:
        f = open(ligandfile, "w")
        f.writelines(l_xyz_lines)
        f.close()
        return ligandfile
    else:
        logging.info("Fatal error: Found multiple ligand for this protein:%s"%proteinfile)
        return False



def get_center(ligand_pdb):
    xyz_lines = open(ligand_pdb,"r").readlines()
    multi_ligand = False
    atom_list = []
    x = y = z = 0
    for xyz_line in xyz_lines:
        if "HETATM" in xyz_line:
            #logging.debug("Check the get center of this protein: %s for this ligand: %s"%(protein_file, ligname))                                                                                                               
            atom_name = xyz_line[12:16]
            if atom_name in atom_list:
                multi_ligand = True
            else:
                atom_list.append(atom_name)
                try:
                    x += float(xyz_line[30:38])
                    y+= float(xyz_line[38:46])
                    z+= float(xyz_line[46:54])
                except:
                    logging.debug("Fatal error: Cannot find the XYZ coordinate for this ligand:%s"%ligand_pdb)
                    return False
    if not multi_ligand:
        lig_center = "%8.3f, %8.3f, %8.3f"%(x/len(atom_list), y/len(atom_list), z/len(atom_list))
        logging.debug("Ligand center for this case:%s is %s"%(ligand_pdb, lig_center))
        return lig_center
    else:
        logging.debug("Fatal error: Found multiple ligands in file:%s"%ligand_pdb)
        return False





        
def main_gendata (s3_result_path, path_2_ent, s4_result_path):
    print 1
    os.chdir(s4_result_path)
    print 5
    current_dir_layer_1 = os.getcwd()
    print 10
    blastnfilterout = glob.glob("%s/*.txt"%s3_result_path)
    summary_file = "%s/summary.txt"%s3_result_path                         
    print 15
    if summary_file in blastnfilterout:                                    
        blastnfilterout.remove(summary_file)
    problem_cases = []
    valid_cases = []
    all_cases = []
    for single_bfout in blastnfilterout:
        print 20
        valid = False
        target_name = os.path.basename(single_bfout).split(".txt")[0]
        all_cases.append(target_name)
        if not os.path.isdir(target_name):
            os.mkdir(target_name)
        os.chdir(target_name)
        commands.getoutput("cp %s ."%single_bfout)
        try:
            info_dic = parse_txt(single_bfout)
        except:
            logging.info("Could not parse this blasternfilter outfile: %s"%(single_bfout))
            continue
        #check if the is a protein start with LMCSS
        if not "LMCSS" in info_dic:
            logging.info("For this query protein: %s, there is no protein sharing the LMCSS ligand with. Not able to generate Docking grid, pass for this case..."%target_name)
            os.chdir (current_dir_layer_1)
            continue
        elif len(info_dic["LMCSS"]) != 2:
            logging.info("For this query protein: %s, the LMCSS protein has wrong number of informations associate with this id..."%target_name)
            os.chdir(current_dir_layer_1)
            continue
        elif not "inchi" in info_dic:
            logging.info("For this query protein: %s, there is no inchi info for the ligand..."%target_name)
            os.chdir(current_dir_layer_1)
            continue
        else:
            logging.info("============Start to work in this query protein: %s============"%target_name)

        #get all ligand related files
        try:
            inchi = info_dic["inchi"]
            ligand_title = info_dic["ligand"]
        #probably will come back sliu
            generate_ligand(inchi, ligand_title)
        except:
            logging.info("Could not generate the ligand for this case: %s"%(target_name))
            os.chdir(current_dir_layer_1)                     
            continue
        #copy all ent file here and rename it
        query_pro = info_dic["query"]
        LMCSS_pro_id = info_dic["LMCSS"][0]
        LMCSS_pro_ligand = info_dic["LMCSS"][1]
        LMCSS_pdb_folder_name = LMCSS_pro_id[1:3]
        LMCSS_ent_file = "pdb" + LMCSS_pro_id  + ".ent"
        LMCSS_pdbloc = os.path.join(path_2_ent, LMCSS_pdb_folder_name, LMCSS_ent_file)
        if not os.path.isfile(LMCSS_pdbloc):            
            logging.info("Unable to find the ent file associate with the LMCSS pdb: %s at location %s"%(LMCSS_pro_id, LMCSS_pdbloc))
            os.chdir(current_dir_layer_1)
            continue
        else:
            #TC = target candidate                                
            TC_id = "%s_%s"%(query_pro, LMCSS_pro_id)             
            LMCSS_protein_name = "LMCSS-%s-%s.pdb"%(TC_id, LMCSS_pro_ligand)
            LMCSS_ligand_name = "LMCSS-%s-%s-lig.pdb"%(TC_id, LMCSS_pro_ligand)
            #try to extract the LMCSS ligand from the LMCSS protein
            commands.getoutput("cp %s %s"%(LMCSS_pdbloc,LMCSS_protein_name))
            if not pull_ligand_out (LMCSS_protein_name, LMCSS_pro_ligand, LMCSS_ligand_name):
                logging.info('Unable to remove ligand from %s. Skipping!' %(LMCSS_protein_name))
                os.chdir(current_dir_layer_1)
                continue
                
            logging.info("Succsessfully generate this protein:%s"%LMCSS_protein_name)

            ## Get the center of mass for the "LMCSS" candidate ligand (all of the other candidates have been aligned to this one)                                                                                          
            #LMCSS_ligand_filename = glob.glob('LMCSS-%s_%s-%s-lig.pdb'%(target_name,LMCSS_pro_id-LMCSS_pro_ligand))
            ligand_center = get_center (LMCSS_ligand_name)
            if not ligand_center:
                logging.info("Unable to find the center of the ligand for the LMCSS candidate pdb: %s"%(pot_target_id))
                os.chdir(current_dir_layer_1)
                continue
            else:
                with open("center.txt" , "w") as center_file:
                    center_file.writelines(ligand_center)





            for rest_protein in ("SMCSS", "hiResHolo", "hiResApo"):
                if rest_protein in info_dic:
                    if len(info_dic[rest_protein]) == 2:
                        rest_protein_id = info_dic[rest_protein][0] 
                        rest_ligand_id = info_dic[rest_protein][1]
                    else:
                        rest_protein_id = info_dic[rest_protein]
                        rest_ligand_id = False
                    rest_protein_folder_name = rest_protein_id[1:3]
                    rest_ent_file = "pdb" + rest_protein_id  + ".ent"       
                    rest_pdbloc = os.path.join(path_2_ent, rest_protein_folder_name, rest_ent_file)
                    if not os.path.isfile(rest_pdbloc):       
                        logging.info("Unable to find the ent file associate with the %s pdb with ID: %s at location %s"%(rest_protein, rest_protein_id, rest_pdbloc))
                        os.chdir(current_dir_layer_1)
                        continue
                    else:
                        TC_id_rest = "%s_%s"%(query_pro, rest_protein_id)
                        if rest_ligand_id:
                            rest_protein_name = "%s-%s-%s.pdb"%(rest_protein, TC_id_rest, rest_ligand_id)
                            commands.getoutput("cp %s %s"%(rest_pdbloc, rest_protein_name))
                        else:
                            rest_protein_name = "%s-%s.pdb"%(rest_protein, TC_id_rest)
                            commands.getoutput("cp %s %s"%(rest_pdbloc, rest_protein_name))
                        #align all rest proteins onto the LMCSS protein
                        
                        #try:
                        do_alignment = align_proteins (LMCSS_protein_name, rest_protein_name, rest_protein_name)
                        if do_alignment == False:
#except:
                            logging.info("The alignment could not be done for this protein:%s"%(rest_protein_name))
                            os.chdir(current_dir_layer_1)
                            continue

                            #remove the "rot-%s"%LMCSS_protein_name
                        if os.path.isfile("rot-%s"%LMCSS_protein_name):
                            commands.getoutput("rm rot-%s"%LMCSS_protein_name) 
                        logging.info("Succsessfully generate this protein:%s"%(rest_protein_name))
                        #here we got the valid aligned structures
                        valid = True
        os.chdir(current_dir_layer_1)        
        #here remove all the folder which don't have full information because of all errors which showed in the final.log
        if valid:
            valid_cases.append(target_name)
    for case in all_cases:
        if case not in valid_cases:
            problem_cases.append(case)
            if not os.path.isdir("error_container"):
                os.mkdir("error_container")
            commands.getoutput("mv %s error_container"%case)
    logging.info("Finish generating the challenge data. The problematic cases are: %s"%problem_cases)

if ("__main__") == (__name__):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-p", "--pdbdb", metavar = "PATH", help = "PDB DATABANK which we will dock into")
    parser.add_option("-c", "--candidatedir", metavar="PATH", help = "PATH where we could find the stage 3 output")
    parser.add_option("-o", "--outdir", metavar = "PATH", help = "PATH where we run stage 4")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.DEBUG )
    (opt, args) = parser.parse_args()
    pdb_location = opt.pdbdb
    stage_3_result = opt.candidatedir
    stage_4_result = opt.outdir
    running_dir = os.getcwd()
    main_gendata(stage_3_result, pdb_location, stage_4_result)
    log_file_path = os.path.join(running_dir, 'final.log')
    commands.getoutput("mv %s %s"%(log_file_path,stage_4_result))
