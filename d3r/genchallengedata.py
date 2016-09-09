#!/usr/bin/env python
__author__ = 'sliu'
#from stage 3 blastnfilter result generating the stage 4 dataset to release to the participant
import os
import sys
import glob
import logging
from rdkit import Chem
from rdkit.Chem import AllChem as allChem
import commands
from Bio import PDB

def split_chain(pdb_filename, out_path, pdb_id, chain_letters):
    parser = PDB.PDBParser()
    writer = PDB.PDBIO()
    chain_letters = [chain.upper() for chain in chain_letters]
    struct = parser.get_structure (pdb_id, pdb_filename)
    writer.set_structure(struct)
    writer.save(out_path, select=SelectChains(chain_letters))

class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

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
        info_value_3 = False
        if line.split(",")[0] == "inchi":
            info_title = line.split(", InChI=")[0]
            info_value_1 = line.split(", InChI=")[1].split("\n")[0]
        elif line.split(",")[0] in ["LMCSS", "SMCSS", "hiResHolo", "hiTanimoto"]:
            info_title = line.split(",")[0]
            info_value_1 = line.split(",")[1].strip()
            info_value_2 = line.split(",")[2].strip()
            info_value_3 = line.split(",")[3].split(":")[1].strip()
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
        if info_value_3:                                                 
            if info_title not in info_dic:                    
                info_dic[info_title] = (info_value_1,info_value_2, info_value_3)           
        else:
            if info_title not in info_dic:
                info_dic[info_title] = (info_value_1,info_value_2) 
    #check if the query protein id is identical with the filename
    if not info_dic["query"] == target_name_txt:
        logging.info("The query pdb ID:%s is different from the filename: %s, need to check"%(info_dic["query"], target_name_txt)) 
    return info_dic  


def align_proteins (target_protein, pre_prepare_protein, post_prepare_protein):
    commands.getoutput("$SCHRODINGER/utilities/structalign %s %s"%(target_protein, pre_prepare_protein))
    rotated_protein = "rot-" + pre_prepare_protein
    if os.path.isfile(rotated_protein):
        commands.getoutput("mv %s %s"%(rotated_protein, post_prepare_protein))
        return True
    else:
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
    #allChem.EmbedMolecule(rd_mol_H)
    #allChem.UFFOptimizeMolecule(rd_mol_H)
    allChem.Compute2DCoords(rd_mol)
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
    if multi_ligand:
        logging.info("Warning: Found multiple ligand for this protein:%s"%proteinfile)
    #here change to save ligand even there are mutiple ligands
    f = open(ligandfile, "w")
    f.writelines(l_xyz_lines)
    f.close()
    return ligandfile

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
    os.chdir(s4_result_path)
    current_dir_layer_1 = os.getcwd()
    blastnfilterout = glob.glob("%s/*.txt"%s3_result_path)
    summary_file = "%s/summary.txt"%s3_result_path                         
    if summary_file in blastnfilterout:                                    
        blastnfilterout.remove(summary_file)
    problem_cases = []
    valid_cases = []
    all_cases = []
    for single_bfout in blastnfilterout:
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
        elif len(info_dic["LMCSS"]) != 3:
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
        LMCSS_mcss_chain = info_dic["LMCSS"][2]
        LMCSS_pdb_folder_name = LMCSS_pro_id[1:3]
        LMCSS_ent_file = "pdb" + LMCSS_pro_id  + ".ent"
        LMCSS_pdbloc = os.path.join(path_2_ent, LMCSS_pdb_folder_name, LMCSS_ent_file)
        if not os.path.isfile(LMCSS_pdbloc):            
            logging.info("Unable to find the ent file associate with the LMCSS pdb: %s at location %s"%(LMCSS_pro_id, LMCSS_pdbloc))
            os.chdir(current_dic_layer_1)
            continue
        else:
            #TC = target candidate                                
            TC_id = "%s_%s"%(query_pro, LMCSS_pro_id)             
            LMCSS_protein_name = "LMCSS-%s-%s.pdb"%(TC_id, LMCSS_pro_ligand)
            LMCSS_ligand_name = "LMCSS-%s-%s-lig.pdb"%(TC_id, LMCSS_pro_ligand)
            #try to extract the LMCSS ligand from the LMCSS protein
            commands.getoutput("cp %s %s"%(LMCSS_pdbloc,LMCSS_protein_name))
            #extract chain where contain the LMCSS mcss
            split_chain(LMCSS_protein_name, LMCSS_protein_name, LMCSS_pro_id, LMCSS_mcss_chain)
            if not pull_ligand_out (LMCSS_protein_name, LMCSS_pro_ligand, LMCSS_ligand_name):
                os.chdir(current_dir_layer_1)
                continue
                
            logging.info("Succsessfully generate this protein:%s"%LMCSS_protein_name)
            ligand_center = get_center (LMCSS_ligand_name)
            if not ligand_center:
                logging.info("Unable to find the center of the ligand for the LMCSS candidate pdb: %s"%(pot_target_id))
                os.chdir(current_dir_layer_1)
                continue
            else:
                with open("center.txt" , "w") as center_file:
                    center_file.writelines(ligand_center)
            for rest_protein in ("SMCSS", "hiResHolo", "hiResApo", "hiTanimoto"):
                if rest_protein in info_dic:
                    if len(info_dic[rest_protein]) == 3:
                        rest_protein_id = info_dic[rest_protein][0] 
                        rest_ligand_id = info_dic[rest_protein][1]
                        rest_ligand_chain = info_dic[rest_protein][2]
                    else:
                        rest_protein_id = info_dic[rest_protein]
                        rest_ligand_id = False
                    rest_protein_folder_name = rest_protein_id[1:3]
                    rest_ent_file = "pdb" + rest_protein_id  + ".ent"       
                    rest_pdbloc = os.path.join(path_2_ent, rest_protein_folder_name, rest_ent_file)
                    if not os.path.isfile(rest_pdbloc):       
                        logging.info("Unable to find the ent file associate with the %s pdb with ID: %s at location %s"%(rest_protein, rest_protein_id, rest_pdbloc))
                        os.chdir(current_dic_layer_1)
                        continue
                    else:
                        TC_id_rest = "%s_%s"%(query_pro, rest_protein_id)
                        if rest_ligand_id:
                            rest_protein_name = "%s-%s-%s.pdb"%(rest_protein, TC_id_rest, rest_ligand_id)
                            commands.getoutput("cp %s %s"%(rest_pdbloc, rest_protein_name))
                            #extract ligand chain from proteins with ligand, original just do this in LMCSS, 0908 sliu 
                            split_chain(rest_protein_name, rest_protein_name, rest_protein_id, rest_ligand_chain)
                        else:
                            rest_protein_name = "%s-%s.pdb"%(rest_protein, TC_id_rest)
                            commands.getoutput("cp %s %s"%(rest_pdbloc, rest_protein_name))
                        #align all rest proteins onto the LMCSS protein
                        
                        try:
                            do_alignment = align_proteins (LMCSS_protein_name, rest_protein_name, rest_protein_name)
                        except:
                            logging.info("The alignment could not be done for this protein:%s"%(rest_protein_name))
                            os.chdir(current_dic_layer_1)
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
