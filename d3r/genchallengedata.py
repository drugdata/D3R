#!/usr/bin/evn python
__author__ = 'sliu'
#from stage 3 blastnfilter result generating the stage 4 dataset to release to the participant
import os
import sys
import glob
import logging
from rdkit import Chem
from rdkit.Chem import AllChem as allChem
import commands


def parse_txt (txt_filename):
    #parsing the result txt file and give back the

    #1, inchi string for the ligand 
    #2, largest, smallest, holo, apo structure ID
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
    allChem.Compute2DCoords(rd_mol)
    Chem.MolToMolFile (rd_mol, mol_filename, includeStereo=True)

def main_gendata (s3_result_path, path_2_ent, s4_result_path):
    os.chdir(s4_result_path)
    current_dir_layer_1 = os.getcwd()
    blastnfilterout = glob.glob("%s/*.txt"%s3_result_path)
    summary_file = "%s/summary.txt"%s3_result_path                         
    if summary_file in blastnfilterout:                                    
        blastnfilterout.remove(summary_file)
    for single_bfout in blastnfilterout:
        target_name = os.path.basename(single_bfout).split(".txt")[0]
        if not os.path.isdir(target_name):
            os.mkdir(target_name)
        os.chdir(target_name)
        commands.getoutput("cp %s ."%single_bfout)
        try:
            info_dic = parse_txt(single_bfout)
        except:
            logging.info("Could not parse this blasternfilter outfile: %s"%(single_bfout))
            continue
        #check if the is a protein start with largest
        if not "largest" in info_dic:
            logging.info("For this query protein: %s, there is no protein sharing the largest ligand with. Not able to generate Docking grid, pass for this case..."%target_name)
            os.chdir (current_dir_layer_1)
            continue
        elif len(info_dic["largest"]) != 2:
            logging.info("For this query protein: %s, the laregest protein has wrong number of informations associate with this id..."%target_name)
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
        #copy all ent file here and rename it
        query_pro = info_dic["query"]
        largest_pro_id = info_dic["largest"][0]
        largest_pro_ligand = info_dic["largest"][1]
        largest_pdb_folder_name = largest_pro_id[1:3]
        largest_ent_file = "pdb" + largest_pro_id  + ".ent"
        largest_pdbloc = os.path.join(path_2_ent, largest_pdb_folder_name, largest_ent_file)
        if not os.path.isfile(largest_pdbloc):            
            logging.info("Unable to find the ent file associate with the largest pdb: %s at location %s"%(largest_pro_id, largest_pdbloc))
            os.chdir(current_dir_layer_1)
            continue
        else:
            #TC = target candidate                                
            TC_id = "%s_%s"%(query_pro, largest_pro_id)             
            largest_protein_name = "largest-%s-%s.pdb"%(TC_id, largest_pro_ligand)
            commands.getoutput("cp %s %s"%(largest_pdbloc,largest_protein_name))
            logging.info("Succsessfully generate this protein:%s"%largest_protein_name)
            for rest_protein in ("smallest", "holo", "apo"):
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
                        #align all rest proteins onto the largest protein
                        
                        try:
                            do_alignment = align_proteins (largest_protein_name, rest_protein_name, rest_protein_name)
                        except:
                            logging.info("The alignment could not be done for this protein:%s"%(rest_protein_name))
                            continue

                            #remove the "rot-%s"%largest_protein_name
                        if os.path.isfile("rot-%s"%largest_protein_name):
                            commands.getoutput("rm rot-%s"%largest_protein_name) 
                        logging.info("Succsessfully generate this protein:%s"%(rest_protein_name))

             
        os.chdir(current_dir_layer_1)        

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
