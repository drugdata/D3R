__author__ = 'robswift'
__project__ = 'blastnfilter'

import sys
import os
import textwrap

def print_usage():
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    try:
        rows, columns = os.popen('stty size', 'r').read().split()
        w.width = int(columns)
    except ValueError:
        pass
    w.initial_indent = ''
    w.subsequent_indent = ''

    usage = "usage: blastnfilter --nonpolymertsv --sequencetsv --crystalpH --pdbblastdb --pdbdb --compinchi --outdir"
    print
    print w.fill(usage)
    print

def interface(args):
    flags = ('--nonpolymertsv', '--sequencetsv', '--crystalpH', '--pdbblastdb', '--pdbdb', '--compinchi', '--outdir',)

    if len(args) == 0:
        print "\nNo arguments specified"
        print_usage()
        sys.exit(1)

    for flag in [x for x in args if '--' in x]:
        if flag not in flags:
            print "\n Unrecognized input flag: %s\n" % flag
            print_usage()
            sys.exit(1)

    parser = ParseArgs(args)
    return parser

class ParseArgs:
    def __init__(self, args):
        self.args = args
        self.required = ('--nonpolymertsv', '--sequencetsv', '--crystalpH', '--pdbblastdb', '--pdbdb', '--compinchi',
                         '--outdir',)

    def get_string(self, input_string):
        if input_string in self.required:
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # the flag wasn't set, so exit
                print "\n %s is required" % input_string
                print_usage()
                sys.exit(1)
            # the flag was set, so check if a value was set, otherwise exit
            try:
                if self.args[index] in self.required:
                    print "\n %s was set but a value was not specified" % input_string
                    print_usage()
                    sys.exit(1)
            except IndexError:
                print "\n %s was set but a value was not specified" % input_string
                print_usage()
                sys.exit(1)
            # return the value that was set
            return format(self.args[index])

class SplitInput:
    def __init__(self, itf):
        self.non_polymer = itf.get_string('--nonpolymertsv')
        self.polymer = itf.get_string('--sequencetsv')
        self.ph = itf.get_string('--crystalpH')
        self.blast_db = itf.get_string('--pdbblastdb')
        self.pdb_path = itf.get_string('--pdbdb')
        self.compinchi = itf.get_string('--compinchi')
        self.out = itf.get_string('--outdir')


