#!/usr/bin/env python

__author__ = 'robswift'

import sys
from d3r.utilities import io_parser
from d3r.utilities import run

def main(argv=[__name__]):
    itf = io_parser.interface(argv[1:])
    options = io_parser.SplitInput(itf)
    run.run(options)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
