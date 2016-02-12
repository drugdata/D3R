__author__ = 'robswift'

import os
import sys
from d3r.utilities.writers import WriteLog
from d3r.utilities.writers import WriteText
from d3r.utilities.analysis import InputAnalysis
from d3r.utilities.analysis import OutputAnalysis

def writer(out_dir, filtrate, log_file=False):
    """
    Writes out a txt file describing queries and their corresponding docking hits and writes out a blastnfilter log
    :param out_dir: (string) the absolute path to the directory where the log file will be written
    :param filtrate: a list of d3r.blast.Query objects
    :param log_file: (boolean) indicates whether or not a log file should be written
    """
    if not os.access(out_dir, os.W_OK):
        print "%s is not writeable" % out_dir
        sys.exit(1)

    txt_writer = WriteText(out_dir)
    log_writer = WriteLog(out_dir)
    if not filtrate.triage:
        txt_writer.write_txt(filtrate)
    if log_file:
        log_writer.write_log(filtrate)
    log_writer.close_file()

def input_analysis(out_dir, queries):
    report = InputAnalysis(queries)
    report.print_to_standard_out()
    report.print_to_file(out_dir)


class OutController(object):
    def __init__(self):
        self.report = OutputAnalysis()

    def set_query(self, query):
        if not query.triage:
            self.report.query = query
            self.report.set_target_dict()

    def print_filter_criteria(self, out_dir):
        self.report.print_filter_criteria(out_dir)

    def print_to_file(self, out_dir):
        self.report.print_to_file(out_dir)