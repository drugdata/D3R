======================================
Drug Design Data Resource CELPP Runner
======================================

.. image:: https://img.shields.io/travis/drugdata/D3R.svg
        :target: https://travis-ci.org/drugdata/D3R.svg?branch=master

.. image:: https://img.shields.io/pypi/v/D3R.svg
        :target: https://pypi.python.org/pypi/D3R


Drug Design Data Resource is a suite of software to enable 
filtering, docking, and scoring of new sequences from wwpdb.


Features
--------

 * Works with Python 2.6, 2.7

Requires
--------

 * argparse
 * lockfile
 * psutil
 * biopython
 * xlsxwriter
 * ftpretty
 * NCBI Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (needed by blastnfilter.py)
 * rdkit (needed by blastnfilter.py and proteinligprep.py)
 * schrodinger (needed by proteinligprep.py and glidedocking.py)

Installation
------------

Pip install is coming, but in the meantime:

.. code:: bash

  git clone https://github.com/drugdata/D3R
  cd D3R
  make install

Usage
-----

Run

.. code:: bash
  
  celpprunner.py --help
