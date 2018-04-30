#!/usr/bin/env python

import argparse
import pickle
import logging
import os
import re
import sys
import math

__author__ = 'sliu'

# pass 1, where contains stage 7 folders, 2,
# contain the stage 8 result, 4 contain genchallenge data final files

RMSD_PICKLE = 'RMSD.pickle'
FINAL_LOG = 'final.log'
CHALL_FINAL_LOG = FINAL_LOG
LMCSS = 'LMCSS'
SMCSS = 'SMCSS'
HI_TANIMOTO = 'hiTanimoto'
HI_RESAPO = 'hiResApo'
HI_RESHOLO = 'hiResHolo'
SUMMARY_TXT = 'summary.txt'
OVERALL_RMSD_PREFIX = 'Overall_RMSD_'
CSV_SUFFIX = '.csv'
NUMBER_BINS = 8
BINSIZE = 1

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
    # TODO Need to find a better way then parsing final.log to get # candidates
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


def get_dock_scores_as_list(pickle_file, ctype=LMCSS):
    """Get list of dock scores for candidate type `ctype`
       from pickle_file `pickle_file`
       The format of the pickle file should be
       { 'ligand': {'ctype': score}}
    :param pickle_file: pickle file to parse
    :param ctype: candidate type
    :returns: list of values, no checking on types could be strings
              ints, floats etc.
    """
    dock_scores = []
    if pickle_file is None:
        logger.error('Pickle file is None')
        return dock_scores

    if ctype is None:
        logger.error('Candidate Type is None')
        return dock_scores

    try:
        p_f = open(pickle_file, "rb")
        p_d = pickle.load(p_f)
        p_f.close()
        for ligand in p_d:
            try:
                value = p_d[ligand][ctype]
                dock_scores.append(value)
            except KeyError:
                logger.error('No ' + ctype +
                             ' for ' + ligand +
                             ' skipping...')
                continue
        return dock_scores
    except Exception:
        logger.exception('Caught Exception just returning empty list')
        return []


def get_list_of_stats(list_of_dock_scores):
    """Extract average for given type
       by getting scores for all ligands with scores
       of type `ctype` in `pickle_file`
       The format of the pickle file should be
       { 'ligand': {'ctype': score}}
    :param list_of_dock_scores: list of dock scores
    :returns: tuple (count, min, max, average, median)
    """
    if list_of_dock_scores is None or len(list_of_dock_scores) is 0:
        logger.error('No dock scores to generate stats from')
        return -1, -1, -1, -1, -1
    try:
        count = len(list_of_dock_scores)

        average = sum(list_of_dock_scores)/count
        if count == 1:
            median = list_of_dock_scores[0]
        else:

            sorted_list = sorted(list_of_dock_scores)
            if count % 2 == 0:
                median_index = count / 2
                median = sum(sorted_list[median_index-1:median_index+1])/2
            else:
                median = sorted_list[int(round(count/2))]
        return (count,
                min(list_of_dock_scores),
                max(list_of_dock_scores),
                average, median)
    except Exception:
        logger.exception('Caught Exception just returning -1')
        return -1, -1, -1, -1, -1


def get_histogram_of_dock_scores(list_of_dock_scores, binsize,
                                 number_bins):
    """creates histogram of dock scores
    :param list_of_dock_scores: list of dock scores
    :param binsize: size of each bin
    :param number_bins: number of bins. If set to say 4 then
                        bin 1 will be values between 0 and < binsize
                        bin 2 will be values between binsize & < binsize x 2
                        bin 3 will be values between binsize x 2 &
                                                     < binsize x 3
                        bin 4 will be values binsize x 3 and greater
    :returns list containing number of scores in each bin
    """
    if number_bins is None or number_bins <= 0:
        logger.error('number_bins is None or <= 0')
        return None
    if binsize <= 0:
        logger.error('binsize is <= 0')
        return None
    if list_of_dock_scores is None or len(list_of_dock_scores) is 0:
        logger.error('No dock scores to generate histogram from')
        return [0] * number_bins

    histo = [0] * number_bins
    for score in list_of_dock_scores:
        thebin = int(math.floor(score/binsize))
        if thebin < 0:
            thebin = 0

        if thebin >= number_bins:
            thebin = number_bins - 1
        histo[thebin] += 1
    return histo


def _get_pickle_paths(path_list):
    """Given a list of evaluation task directories
       find all pickle files and return that as
       a list as well as a list of those that failed
    """
    pickle_list = []
    no_pickle_list = []

    if path_list is None:
        logger.error('Path list is None')
        return pickle_list, no_pickle_list

    logger.info('Given ' + str(len(path_list)) + ' paths to look for'
                                                 ' pickle files in')

    for all_stage_7 in path_list:
        full_path = os.path.join(all_stage_7, RMSD_PICKLE)
        if os.path.isfile(full_path):
            pickle_list.append(full_path)
        else:
            logging.info("The pickle file " + full_path +
                         "does not exist")
            no_pickle_list.append(full_path)

    return pickle_list, no_pickle_list


def _get_submission_name_from_pickle_path(path, prefix, suffix,
                                          max_submission_name_width=None):
    """Given path to pickle path in evaluation task
       get the submission name which is basically
       the prefix and suffix values lopped off
    """
    dir_name = os.path.dirname(path)
    base_name = os.path.basename(dir_name)
    sub_name_suffix = re.sub('^' + prefix, '', base_name)
    res = re.sub(suffix, '', sub_name_suffix)
    if max_submission_name_width is None or max_submission_name_width <= 3:
        return res
    if len(res) > max_submission_name_width:
        return res[0:max_submission_name_width-2] + '..'
    return res


def generate_overall_csv(evaluation_path, challenge_dir, post_evaluation_path,
                         candidates_type=LMCSS,
                         eval_stage_prefix='',
                         eval_suffix='',
                         submission_name_width=30):

    not_valid_pickle = 0

    all_pickle_files, non_pickle_list = _get_pickle_paths(evaluation_path)

    candidates_report = os.path.join(challenge_dir, CHALL_FINAL_LOG)

    # Doing 2 checks cause older versions had mispelled Succsessfully in
    # final.log
    # output and new ones fixed this.
    mispelled_candidates = check_case_number(candidates_report,
                                             "Succsessfully generate " +
                                             "this protein:" + candidates_type)
    correctspelled_candidates = check_case_number(candidates_report,
                                                  "Successfully generate "
                                                  "this protein:" +
                                                  candidates_type)
    total_candidates = 0
    if mispelled_candidates > 0:
        total_candidates += mispelled_candidates
    if correctspelled_candidates > 0:
        total_candidates += correctspelled_candidates

    overall_csv = open(os.path.join(post_evaluation_path,
                                    OVERALL_RMSD_PREFIX + candidates_type +
                                    CSV_SUFFIX), 'w')

    summary_txt = open(os.path.join(post_evaluation_path,
                                    SUMMARY_TXT), 'a')

    s_line_format = ('%-' + str(submission_name_width+1) +
                     's%-10s%-10s%-10s%-10s%-12s'
                     '%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s\n')
    headerline = (s_line_format % (candidates_type + ' (' +
                                   str(total_candidates) +
                                   ' dockable)',
                                   '#Docked',
                                   'Min RMSD',
                                   'Max RMSD',
                                   'Mean RMSD',
                                   'Median RMSD',
                                   '0<1',
                                   '1<2',
                                   '2<3',
                                   '3<4',
                                   '4<5',
                                   '5<6',
                                   '6<7',
                                   '7+'))
    summary_txt.write(headerline)

    summary_txt.write('-'*int(len(headerline)-1) + '\n')

    data_line_format = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'
    full_data_lines = [data_line_format % ("SubmissionID for " +
                                           candidates_type,
                                           "# Docked",
                                           "# Dockable",
                                           "Min RMSD",
                                           "Max RMSD",
                                           "Mean RMSD",
                                           "Median RMSD",
                                           "0<1",
                                           "1<2",
                                           "2<3",
                                           "3<4",
                                           "4<5",
                                           "5<6",
                                           "6<7",
                                           "7+")]
    for p_file in all_pickle_files:
        sub_name = _get_submission_name_from_pickle_path(p_file,
                                                         eval_stage_prefix,
                                                         eval_suffix)
        tsub_name = _get_submission_name_from_pickle_path(p_file,
                                                          eval_stage_prefix,
                                                          eval_suffix,
                                                          max_submission_name_width=submission_name_width)
        try:
            d_scores = get_dock_scores_as_list(p_file,
                                               ctype=candidates_type)
            (num_docked, minrmsd, maxrmsd,
             avgrmsd, medianrmsd) = get_list_of_stats(d_scores)

            histo = get_histogram_of_dock_scores(d_scores, BINSIZE,
                                                 NUMBER_BINS)
            if histo is None:
                histo = [0] * NUMBER_BINS

            pc_docked = '(NA%)'
            if total_candidates > 0:
                percent_docked = float(num_docked) / float(total_candidates)
                percent_docked *= 100.0
                pc_docked = ' (%3.0f%%)' % percent_docked

            full_data_lines.append(data_line_format % (sub_name,
                                                       num_docked,
                                                       total_candidates,
                                                       '%.2f' % minrmsd,
                                                       '%.2f' % maxrmsd,
                                                       '%.2f' % avgrmsd,
                                                       '%.2f' % medianrmsd,
                                                       '%d' % histo[0],
                                                       '%d' % histo[1],
                                                       '%d' % histo[2],
                                                       '%d' % histo[3],
                                                       '%d' % histo[4],
                                                       '%d' % histo[5],
                                                       '%d' % histo[6],
                                                       '%d' % histo[7]))
            summary_txt.write(s_line_format % (tsub_name,
                                               str(num_docked) + pc_docked,
                                               '%.2f' % minrmsd,
                                               '%.2f' % maxrmsd,
                                               '%.2f' % avgrmsd,
                                               '%.2f' % medianrmsd,
                                               '%d' % histo[0],
                                               '%d' % histo[1],
                                               '%d' % histo[2],
                                               '%d' % histo[3],
                                               '%d' % histo[4],
                                               '%d' % histo[5],
                                               '%d' % histo[6],
                                               '%d' % histo[7]))
        except IOError:
            logging.exception('Error writing entry')
        except (ValueError, TypeError):
            full_data_lines.append(data_line_format % (sub_name,
                                                       "N/A",
                                                       total_candidates,
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A",
                                                       "N/A"))
            summary_txt.write(s_line_format % (tsub_name,
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A",
                                               "N/A"))
            logging.exception("This pickle file :%s do not have valid "
                              "data for "
                              "candidates type  %s" % (p_file,
                                                       candidates_type))
            not_valid_pickle += 1

    for np_file in non_pickle_list:
        sub_name = _get_submission_name_from_pickle_path(np_file,
                                                         eval_stage_prefix,
                                                         eval_suffix)
        tsub_name = _get_submission_name_from_pickle_path(np_file,
                                                          eval_stage_prefix,
                                                          eval_suffix,
                                                          max_submission_name_width=submission_name_width)
        full_data_lines.append(data_line_format % (sub_name,
                                                   "N/A",
                                                   total_candidates,
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A",
                                                   "N/A"))
        summary_txt.write(s_line_format % (tsub_name,
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A",
                                           "N/A"))

    logging.info("We got : %s cases without pickle files and %s cases have "
                 "invalid pickle files for "
                 "candidates type %s" % (str(len(non_pickle_list)),
                                         not_valid_pickle,
                                         candidates_type))
    overall_csv.writelines(full_data_lines)
    overall_csv.close()

    summary_txt.write('\n\n')
    summary_txt.close()


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
    parser.add_argument("--histogrambinsize", default=1,
                        help='Sets histogram bin size in angstroms '
                             '(default 1)')
    parser.add_argument("--histogrambincount", default=3)
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
    for candidate_type in [LMCSS, SMCSS, HI_TANIMOTO,
                           HI_RESAPO, HI_RESHOLO]:
        try:
            generate_overall_csv(opt.evaluationdir,
                                 opt.challengedir,
                                 opt.outdir,
                                 candidates_type=candidate_type,
                                 eval_stage_prefix=opt.stageprefix,
                                 eval_suffix=opt.evaluationsuffix)
        except Exception:
            logger.exception('Caught exception running generate_overall_csv '
                             'for candidate type ' + candidate_type)
            continue
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
