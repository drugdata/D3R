__author__ = 'robswift'

import os
import sys
from writers import WriteLog
from writers import WriteTxt
from d3r.utilities.analysis import InputAnalysis
from d3r.utilities.analysis import OutputAnalysis

def writer(out_dir, filtrate, log_file = False):
    """
    Writes out a txt file describing queries and their corresponding docking hits and writes out a blastnfilter log
    :param out_dir: (string) the absolute path to the directory where the log file will be written
    :param filtrate: a list of d3r.blast.Query objects
    """
    if not os.access(out_dir, os.W_OK):
        print "%s is not writeable" % out_dir
        sys.exit(1)

    txt_writer = WriteTxt(out_dir)
    log_writer = WriteLog(out_dir)
    if not filtrate.triage:
        txt_writer.write_txt(filtrate)
    if log_file:
        log_writer.write_log(filtrate)
    log_writer.close_file()


#    print "Target ID: %s" % filtrate.pdb_id
#    if filtrate.triage:
#        for reason in filtrate.reasons_to_triage:
#            print "Triage Reason: %s" % reason
#    for hit in filtrate.hits:
#        print "  Test ID: %s" % hit.pdb_id
#        if hit.triage:
#            for reason in hit.reasons_to_triage:
#                print "  Triage Reason: %s" % reason
#        if hit.retain:
#            for reason in hit.reasons_to_retain:
#                print "  Retain Reason: %s" % reason
#        print "  Method: %s" % hit.exp_method
#        print "  Resolution: %s" % hit.resolution
#        print "  Chain count: %s" % hit.chain_count
#        print "  Dock count: %s" % hit.dock_count
#        for seq in hit.sequences:
#            for qa in seq.query_alignments:
#                print "  Chain_%s  Identity: %.2f (chain_%s)|Coverage: %.2f (chain_%s)" % (seq.hit_chain_id,
#                                                                                           qa.identity,
#                                                                                           qa.query_chain_id,
#                                                                                           qa.coverage,
#                                                                                           qa.query_chain_id)
#        print ""

def input_analysis(out_dir, queries):
    report = InputAnalysis(queries)
    report.print_to_standard_out()
    report.print_to_file(out_dir)


class OutController(object):
    def __init__(self):
        self.report = OutputAnalysis()

    def set_target(self, query):
        if not query.triage:
            self.report.query = query
            self.report.set_target_dict()

    def print_filter_criteria(self, out_dir):
        self.report.print_filter_criteria(out_dir)

    def print_to_file(self, out_dir):
        self.report.print_to_file(out_dir)