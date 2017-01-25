#!/usr/bin/env python

__author__ = 'sliu'
import pickle
import logging
import commands
import glob
import os
import sys
#pass 1, where contains stage 7 folders, 2, contain the stage 8 result, 4 contain genchallenge data final files

def check_case_number (log_file, phase):
    log_f = open(log_file, "r")
    log_lines = log_f.readlines()
    total_number = 0
    for log_line in log_lines:
        if phase in log_line:
            total_number += 1
    return total_number
def extract_ave (pickle_file, candidate_type = "LMCSS"):
    p_f = open(pickle_file, "r")
    p_d = pickle.load(p_f)
    p_f.close()
    data = []
    for ligand in p_d:
        try:
            value = p_d[ligand][candidate_type]
            data.append(value)
        except:
            continue
    number_of_bins = len(data)
    average = sum(data)/number_of_bins
    return number_of_bins, average

def generate_overall_csv (evaluation_path, challenge_dir, post_evaluation_path, candidates_type = "LMCSS"):
    #evaluation_path =["stage.7.glide.evaluation", "stage.7.autodockvina.evaluation" ... ]
    all_pickle_files = []
    non_pickle_case = 0
    not_valid_pickle = 0
    for all_stage_7 in evaluation_path:
        if os.path.isfile("%s/RMSD.pickle"%all_stage_7):
            all_pickle_files.append("%s/RMSD.pickle"%all_stage_7)
        else:
            logging.info("The pickle file :%s/RMSD.pickle is not exist"%all_stage_7)
            non_pickle_case += 1
    #all_pickle_files = glob.glob("%s/stage.7.*/RMSD.pickle"%evaluation_path)
    overall_csv = open("%s/Overall_RMSD_%s.csv"%(post_evaluation_path, candidates_type), "w")
    candidates_report = os.path.join(challenge_dir, "final.log")
    total_candidates = check_case_number(candidates_report, "Succsessfully generate this protein:%s"%candidates_type)
    full_data_lines = ["%-30s,%-30s,%-30s,%-30s\n"%("Submission ID", "Number docked", "Number to be docked", "Ave RMSD for %s"%(candidates_type))]
    for submitted_pickle in all_pickle_files:
        submission_name = os.path.splitext(os.path.dirname(submitted_pickle))[0].split("stage.7.")[1]
        try:
            number_of_bins, ave = extract_ave(submitted_pickle, candidate_type = candidates_type)
            full_data_lines.append("%-30s,%-30s,%-30s,%-30s\n"%(submission_name, number_of_bins, total_candidates, ave))
        except:
            full_data_lines.append("%-30s,%-30s,%-30s,%-30s\n"%(submission_name, "N/A", total_candidates, "N/A"))
            logging.info("This pickle file :%s do not have valid data for candidates type %s"%(submitted_pickle, candidates_type))
            not_valid_pickle += 1
    logging.info("We got : %s cases without pickle files and %s cases have invalid pickle files for candidates type %s"%(non_pickle_case, not_valid_pickle, candidates_type))
    overall_csv.writelines(full_data_lines)
    overall_csv.close()

if ("__main__") == (__name__):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-e", "--evaluationdir", metavar="PATH", action='append',
                      help="PATH where we could find the evaluation stage output")

    parser.add_argument("-o", "--outdir", metavar="PATH",
                      help="PATH where we will run the evaluate stage")

    parser.add_argument("-c", "--chanllengedir", metavar="PATH",
                      help="PATH to celpp week directory under challenge task ")       
    opt = parser.parse_args()
    challengeDir = opt.chanllengedir
    postDir = opt.outdir
    evaluateDir = opt.evaluationdir

    logger = logging.getLogger()
    logging.basicConfig(format='%(asctime)s: %(message)s',
                        datefmt='%m/%d/%y %I:%M:%S', filename='%s/final.log'%postDir,
                        filemode='w', level=logging.INFO)
    #main calculation
    for candidate_type in ["LMCSS", "SMCSS", "hiTanimoto", "hiResApo", "hiResHolo"]:
        try:
            generate_overall_csv(evaluateDir, challengeDir, postDir, candidates_type = candidate_type)
            #commands.getoutput("mv Overall_RMSD_%s.csv %s"%(candidate_type,postDir) )
        except:
            continue

