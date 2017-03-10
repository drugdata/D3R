#!/usr/bin/env python

import argparse
import pickle
import logging
import os
import re
import sys

__author__ = 'sliu'

# pass 1, where contains stage 7 folders, 2,
# contain the stage 8 result, 4 contain genchallenge data final files

RMSD_PICKLE = 'RMSD.pickle'
FINAL_LOG = 'final.log'
CHALL_FINAL_LOG = FINAL_LOG
LMCSS = 'LMCSS'

logger = logging.getLogger()


def check_case_number(log_file, phrase):
    """Counts number of occurrences of `phrase`
       in file path passed in.
    :param log_file: path to file to examine
    :param phrase: phrase to search for in each line of file
    :returns: count of occurrences as int upon success, -1 if
              there was an error or if `log_file` or 1phrase`
              is None
    """

    total_number = -1
    if log_file is None:
        logger.error('log file is None')
        return total_number

    if phrase is None:
        logger.error('phrase is None')
        return total_number

    log_f = None
    try:
        log_f = open(log_file, "r")
        total_number = 0
        for log_line in log_f:
            if phrase in log_line:
                total_number += 1
    except IOError:
        logger.exception('Caught exception attempting to load file: ' +
                         log_file)
    finally:
        if log_f is not None:
            log_f.close()
    return total_number


def get_numdocked_and_average(pickle_file, ctype=LMCSS):
    """Extract average for given type
       by getting scores for all ligands with scores
       of type `ctype` in `pickle_file`
       The format of the pickle file should be
       { 'ligand': {'ctype': score}}
    :param pickle_file: pickle file to parse
    :param ctype: candidate type
    :returns: tuple (# docked, average) unless there is an
              error with arguments or reading pickle file then
              (-1, -1) is returned. If no docked ligands are
              found then (0, -1) is returned
    """
    if pickle_file is None:
        logger.error('Pickle file is None')
        return -1, -1

    if ctype is None:
        logger.error('Candidate Type is None')
        return -1, -1

    try:
        p_f = open(pickle_file, "rb")
        p_d = pickle.load(p_f)
        p_f.close()
        data = []
        for ligand in p_d:
            try:
                value = p_d[ligand][ctype]
                data.append(value)
            except KeyError:
                logger.error('No ' + ctype +
                             ' for ' + ligand +
                             ' skipping...')
                continue
        number_of_bins = len(data)

        if number_of_bins <= 0:
            logger.error('Zero docked results for'
                         ' pickle file: ' + pickle_file +
                         ' returning 0, -1 cause it is'
                         ' not possible to calculate average')
            return 0, -1

        average = sum(data)/number_of_bins
        return number_of_bins, average
    except Exception:
        logger.exception('Caught Exception just returning -1, -1')
        return -1, -1


def _get_pickle_paths_from_evaluation_directories(path_list):
    """Given a list of evaluation task directories
       find all pickle files and return that as
       a list as well as a count of those that failed
    """
    pickle_list = []
    no_pickle_count = 0
    for all_stage_7 in path_list:
        full_path = os.path.join(all_stage_7, RMSD_PICKLE)
        if os.path.isfile(full_path):
            pickle_list.append(full_path)
        else:
            logging.info("The pickle file " + full_path +
                         "does not exist")
            no_pickle_count += 1

    return pickle_list, no_pickle_count


def generate_overall_csv(evaluation_path, challenge_dir, post_evaluation_path,
                         candidates_type=LMCSS,
                         eval_stage_prefix='stage.7.',
                         eval_suffix='.evaluation'):

    not_valid_pickle = 0

    all_pickle_files, non_pickle_case = _get_pickle_paths_from_evaluation_directories(evaluation_path)

    candidates_report = os.path.join(challenge_dir, CHALL_FINAL_LOG)
    total_candidates = check_case_number(candidates_report,
                                         "Succsessfully generate "
                                         "this protein:" + candidates_type)

    overall_csv = open(os.path.join(post_evaluation_path,
                                    'Overall_RMSD_' + candidates_type +
                                    '.csv'),
                       'w')
    # TODO Add generation of .txt file that we can put
    # TODO into summary email that is sent out
    data_line_format = '%-30s,%-10s,%-12s,%-22s\n'
    full_data_lines = [data_line_format % ("Submission ID",
                                                      "# Docked",
                                                      "# Dockable",
                                                      "Ave RMSD "
                                                      "for %s" %
                                                      candidates_type)]
    for p_file in all_pickle_files:
        dir_name = os.path.dirname(p_file)
        base_name = os.path.basename(dir_name)
        sub_name_suffix = re.sub('^' + eval_stage_prefix, '', base_name)
        sub_name = re.sub(eval_suffix, '', sub_name_suffix)
        try:
            num_docked, avg = get_numdocked_and_average(p_file,
                                                        ctype=candidates_type)
            full_data_lines.append(data_line_format % (sub_name,
                                                       num_docked,
                                                       total_candidates, avg))
        except:
            full_data_lines.append(data_line_format % (sub_name,
                                                       "N/A",
                                                       total_candidates,
                                                       "N/A"))
            logging.info("This pickle file :%s do not have valid data for "
                         "candidates type  %s" % (p_file,
                                                 candidates_type))
            not_valid_pickle += 1
    logging.info("We got : %s cases without pickle files and %s cases have "
                 "invalid pickle files for "
                 "candidates type %s" % (non_pickle_case,
                                         not_valid_pickle,
                                         candidates_type))
    overall_csv.writelines(full_data_lines)
    overall_csv.close()


def main(args):
    """Main entry point of script"""
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

    opt = parser.parse_args(args[1:])
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
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
