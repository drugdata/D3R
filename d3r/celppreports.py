#! /usr/bin/env python

import sys
import os
import argparse
import logging

import d3r
from d3r.celpp import util
from d3r.celpp.task import D3RParameters
from d3r.celpp.blastnfilter import BlastNFilterTask

# create logger
logger = logging.getLogger('d3r.celppreports')
LOG_FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"


def _setup_logging(theargs):
    """Sets up the logging for application

    Loggers are setup for:
    d3r.celppreports
    d3r.celpp.blastnfilter
    d3r.celpp.dataimport
    d3r.celpp.glide
    d3r.celpp.makeblastdb
    d3r.celpp.proteinligprep
    d3r.celpp.evaluation
    d3r.celpp.task
    d3r.celpp.util
    d3r.celpp.uploader

    NOTE:  If new modules are added please add their loggers to this
    function

    The loglevel is set by theargs.loglevel and the format is set
    by LOG_FORMAT set at the top of this module.
    :param: theargs should have .loglevel set to one of the
    following strings: DEBUG, INFO, WARNING, ERROR, CRITICAL
    """
    theargs.logformat = LOG_FORMAT
    theargs.numericloglevel = logging.NOTSET
    if theargs.loglevel == 'DEBUG':
        theargs.numericloglevel = logging.DEBUG
    if theargs.loglevel == 'INFO':
        theargs.numericloglevel = logging.INFO
    if theargs.loglevel == 'WARNING':
        theargs.numericloglevel = logging.WARNING
    if theargs.loglevel == 'ERROR':
        theargs.numericloglevel = logging.ERROR
    if theargs.loglevel == 'CRITICAL':
        theargs.numericloglevel = logging.CRITICAL

    logger.setLevel(theargs.numericloglevel)
    logging.basicConfig(format=theargs.logformat)

    # There should be a line below for every package aka .py file
    # under celpp module
    logging.getLogger('d3r.celpp.blastnfilter')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.dataimport').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.glide').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.makeblastdb')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.proteinligprep')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.evaluation').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.task').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.util').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.uploader').setLevel(theargs.numericloglevel)


def generate_reports(theargs):
    """Generates reports
    """
    celpp_years = util.get_all_celpp_years(theargs.celppdir)

    if theargs.outdir is None:
        raise Exception('--outdir must be set')

    if not os.path.isdir(theargs.outdir):
        os.makedirs(theargs.outdir)

    f = open(os.path.join(theargs.outdir, 'blastnfilter.summary.csv'), 'w')
    f.write('Week #, Year, Complexes, Dockable complexes, Dockable monomers, '
            'Targets Found\n')
    for year in celpp_years:
        logger.info('Examining year ' + year)
        for week in util.get_all_celpp_weeks(os.path.join(theargs.celppdir,
                                                          year)):
            logger.debug('Examining week ' + week)
            the_dir = os.path.join(theargs.celppdir, year,
                                   util.DATA_SET_WEEK_PREFIX + week)
            blast = BlastNFilterTask(the_dir, theargs)
            summary = blast.get_blastnfilter_summary()
            f.write(summary.get_csv() + '\n')
    f.flush()
    f.close()


def _parse_arguments(desc, args):
    """Parses command line arguments using argparse.
    """
    parsed_arguments = D3RParameters()

    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("--outdir",
                        help='Directory to write output reports to ')
    parser.add_argument("celppdir",
                        help='Directory where celpp yearly runs reside')
    parser.add_argument("--log", dest="loglevel", choices=['DEBUG',
                        'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level",
                        default='WARNING')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + d3r.__version__))
    return parser.parse_args(args, namespace=parsed_arguments)


def main():
    desc = """
              Generates reports on CELPP weekly runs
              """

    theargs = _parse_arguments(desc, sys.argv[1:])
    theargs.program = sys.argv[0]
    theargs.version = d3r.__version__

    _setup_logging(theargs)

    try:
        generate_reports(theargs)
    except Exception:
        logger.exception("Error caught exception")
        sys.exit(2)

if __name__ == '__main__':
    main()
