#!/usr/bin/env python

__author__ = 'robswift'

import sys
from d3r.utilities import io_parser
from d3r.utilities import run
import logging

# create logger
logger = logging.getLogger('d3r.blastnfilter')
LOG_FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"


def _setup_logging(loglevel, logformat):
    """Sets up logging for the application

    Loggers are setup for
    d3r.blastnfilter
    d3r.blast.ligand

    :param: loglevel logging level which should be set to a valid input
            for logging.setLevel()
    :param: logformat format to set for logging messages
    """
    logger.setLevel(loglevel)
    logging.basicConfig(format=logformat)
    logging.getLogger('d3r.blast.ligand').setLevel(loglevel)
    logging.getLogger('d3r.utilities.run').setLevel(loglevel)


def main(argv=[__name__]):

    itf = io_parser.interface(argv[1:])
    options = io_parser.SplitInput(itf)
    _setup_logging(logging.DEBUG, LOG_FORMAT)
    logger.debug('Starting run ')
    run.run(options)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
