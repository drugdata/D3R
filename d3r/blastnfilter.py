#!/usr/bin/env python

__author__ = 'robswift'

import sys
import logging
import argparse

from d3r.utilities import run
import d3r
from d3r.celpp import util

# create logger
logger = logging.getLogger('d3r.blastnfilter')
LOG_FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"

class CommandLineParameters(object):
    """Holds command line arguments
    """
    pass


def _parse_arguments(desc, args):
    """Parses command line arguments using argparse.
    """
    parsed_arguments = CommandLineParameters()

    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("--outdir", required=True, dest="out",
                        help='Directory to write output to ')
    parser.add_argument("--nonpolymertsv", required=True,
                        dest="non_polymer",
                        help='Nonpolymer tsv file ')
    parser.add_argument("--sequencetsv", dest="polymer",
                        required=True,
                        help='Sequence tsv file ')
    parser.add_argument("--crystalpH", dest="ph",
                        required=True,
                        help='Crystal ph file')
    parser.add_argument("--pdbblastdb", dest="blast_db",
                        required=True,
                        help='Blast database directory')
    parser.add_argument("--pdbdb", dest="pdb_path",
                        required=True,
                        help='PDB Database directory')
    parser.add_argument("--compinchi", required=True,
                        help='Components.inchi file')
    parser.add_argument("--log", dest="loglevel", choices=['DEBUG',
                        'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level",
                        default='WARNING')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + d3r.__version__))
    return parser.parse_args(args, namespace=parsed_arguments)

def main(args):
    desc = """
              Blastnfilter
              """
    options = _parse_arguments(desc, args[1:])
    options.program = args[0]
    options.version = d3r.__version__

    util.setup_logging(options)
    logger.debug('Starting run ')
    try:
        run.run(options)
    except Exception:
        logger.exception("Error caught exception")
        sys.exit(2)
    sys.exit(0)

if __name__ == "__main__":
    main(sys.argv)
