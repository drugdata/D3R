#!/usr/bin/env python

import argparse
import pickle
import logging
import os
import re

__author__ = 'sliu'

# pass 1, where contains stage 7 folders, 2,
# contain the stage 8 result, 4 contain genchallenge data final files

RMSD_PICKLE = 'RMSD.pickle'
FINAL_LOG = 'final.log'
CHALL_FINAL_LOG = FINAL_LOG
LMCSS = 'LMCSS'


def check_case_number(log_file, phase):
    total_number = -1
    log_f = None
    try:
        log_f = open(log_file, "r")
        total_number = 0
        for log_line in log_f:
            if phase in log_line:
                total_number += 1
    except IOError:
        logger.exception('Caught exception attempting to load file: ' +
                         log_file)
    finally:
        if log_f is not None:
            log_f.close()
    return total_number


def extract_ave(pickle_file, ctype=LMCSS):
    """Extract average for given type
    :param pickle_file: pickle file to parse
    :param ctype: candidate type
    :returns: tuple (number of binds, average)
    """
    p_f = open(pickle_file, "r")
    p_d = pickle.load(p_f)
    p_f.close()
    data = []
    for ligand in p_d:
        try:
            value = p_d[ligand][ctype]
            data.append(value)
        except:
            continue
    number_of_bins = len(data)
    average = sum(data)/number_of_bins
    return number_of_bins, average


def generate_overall_csv(evaluation_path, challenge_dir, post_evaluation_path,
                         candidates_type=LMCSS,
                         eval_stage_prefix='stage.7.',
                         eval_suffix='.evaluation'):
    all_pickle_files = []
    non_pickle_case = 0
    not_valid_pickle = 0
    for all_stage_7 in evaluation_path:
        full_path = os.path.join(all_stage_7, RMSD_PICKLE)
        if os.path.isfile(full_path):
            all_pickle_files.append(full_path)
        else:
            logging.info("The pickle file " + full_path +
                         "does not exist")
            non_pickle_case += 1
    overall_csv = open(os.path.join(post_evaluation_path,
                                    'Overall_RMSD_' + candidates_type +
                                    '.csv'),
                       'w')
    candidates_report = os.path.join(challenge_dir, CHALL_FINAL_LOG)
    total_candidates = check_case_number(candidates_report,
                                         "Succsessfully generate "
                                         "this protein:" + candidates_type)

    data_line_format = '%-30s,%-10s,%-12s,%-22s\n'
    full_data_lines = [data_line_format % ("Submission ID",
                                                      "# Docked",
                                                      "# Dockable",
                                                      "Ave RMSD "
                                                      "for %s" %
                                                      candidates_type)]
    for submitted_pickle in all_pickle_files:
        dir_name = os.path.dirname(submitted_pickle)
        base_name = os.path.basename(dir_name)
        sub_name_suffix = re.sub('^' + eval_stage_prefix, '', base_name)
        sub_name = re.sub(eval_suffix, '', sub_name_suffix)
        try:
            number_of_bins, ave = extract_ave(submitted_pickle,
                                              ctype=candidates_type)
            full_data_lines.append(data_line_format % (sub_name,
                                                       number_of_bins,
                                                       total_candidates, ave))
        except:
            full_data_lines.append(data_line_format % (sub_name,
                                                       "N/A",
                                                       total_candidates,
                                                       "N/A"))
            logging.info("This pickle file :%s do not have valid data for "
                         "candidates type %s" % (submitted_pickle,
                                                 candidates_type))
            not_valid_pickle += 1
    logging.info("We got : %s cases without pickle files and %s cases have "
                 "invalid pickle files for "
                 "candidates type %s" % (non_pickle_case,
                                         not_valid_pickle,
                                         candidates_type))
    overall_csv.writelines(full_data_lines)
    overall_csv.close()


if "__main__" == __name__:

    desc = """
              Summarizes results for all docking evaluations
              passed into script
    """
    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("outdir",
                        help="PATH where we will run the evaluate stage")
    parser.add_argument("-e", "--evaluationdir", metavar="PATH",
                        action='append',
                        help="PATH where we could find the evaluation "
                             "stage output")
    parser.add_argument("-c", "--challengedir", metavar="PATH",
                        help="PATH to celpp week directory under challenge "
                             "task ")
    parser.add_argument("--stageprefix", default='stage.7.',
                        help="Evaluation task dir prefix (default stage.7.) ")
    parser.add_argument("--evaluationsuffix",
                        default='.extsubmission.evaluation$|.evaluation$',
                        help="Evaluation task dir suffix regex"
                             "(default .extsubmission.evaluati"
                             "on$|.evaluation$)")

    opt = parser.parse_args()
    logger = logging.getLogger()
    logging.basicConfig(format='%(asctime)s: %(message)s',
                        datefmt='%m/%d/%y %I:%M:%S',
                        filename=os.path.join(opt.outdir, FINAL_LOG),
                        filemode='w', level=logging.INFO)
    # main calculation
    for candidate_type in [LMCSS, "SMCSS", "hiTanimoto",
                           "hiResApo", "hiResHolo"]:
        try:
            generate_overall_csv(opt.evaluationdir,
                                 opt.challengedir,
                                 opt.outdir,
                                 candidates_type=candidate_type,
                                 eval_stage_prefix=opt.stageprefix,
                                 eval_suffix=opt.evaluationsuffix)
        except:
            logger.exception('Caught exception running generate_overall_csv '
                             'for candidate type ' + candidate_type)
            continue
