===============================
Drug Design Data Resource CELPP Runner
===============================

.. image:: https://img.shields.io/travis/drugdata/D3R.svg
        :target: https://travis-ci.org/drugdata/D3R

.. image:: https://img.shields.io/pypi/v/D3R.svg
        :target: https://pypi.python.org/pypi/D3R


Drug Design Data Resource is a suite of software to enable 
filtering, docking, and scoring of new sequences from wwpdb.


Features
--------

 * Works with Python 2.6, 2.7, 3.0, 3.2, 3.4

Requires
--------

 * argparse
 * lockfile
 * psutil
 * biopython
 * xlsxwriter
 * NCBI Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
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
  
  usage: celpprunner.py [-h] [--blastdir BLASTDIR] [--email EMAIL]
                      [--createweekdir] [--customweekdir] --stage STAGE
                      [--blastnfilter BLASTNFILTER]
                      [--postanalysis POSTANALYSIS]
                      [--proteinligprep PROTEINLIGPREP] [--glide GLIDE]
                      [--evaluation EVALUATION] [--pdbdb PDBDB]
                      [--compinchi COMPINCHI] [--pdbfileurl PDBFILEURL]
                      [--log {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                      [--smtp SMTP] [--smtpport SMTPPORT] [--version]
                      celppdir
  
              Runs the 5 stages (import, blast, proteinligprep, glide,
              & evaluation) of CELPP processing pipeline
              (http://www.drugdesigndata.org)
  
              CELPP processing pipeline relies on a set of directories
              with specific structure. The pipeline runs a set of stages
              Each stage has a numerical value and a name. The numerical
              value denotes order and the stage name identifies separate
              tasks to run in the stage.
  
              The filesystem structure of the stage is:
  
              stage.<stage number>.<task name>
  
              The stage(s) run are defined via the required --stage flag.
  
              To run multiple stages serially just pass a comma delimited
              list to the --stage flag. Example: --stage import,blast
  
              NOTE:  When running multiple stages serially the program will
                     not run subsequent stages if a task in a stage fails.
                     Also note order matters, ie putting blast,import will
                     cause celpprunner.py to run blast stage first.
  
              This program drops a pid lockfile
              (celpprunner.<stage>.lockpid) in celppdir to prevent duplicate
              invocation.
  
              When run, this program will examine the stage and see
              if work can be done.  If stage is complete or previous
              steps have not completed, the program will exit silently.
              If previous steps have failed or current stage already
              exists in an error or uncomplete state then program will
              report the error via email using addresses set in --email
              flag. Errors will also be reported via stderr/stdout.
              The program will also exit with nonzero exit code.
  
              This program utilizes simple token files to denote stage
              completion.  If within the stage directory there is a:
  
              'complete' file - then stage is done and no other
                                checking is done.
  
              'error' file - then stage failed.
  
              'start' file - then stage is running.
  
              Notification of stage start and end will be sent to
              addresses set via --email flag.
  
              Unless --customweekdir is set, this program will
              examine the 'celppdir' (last argument passed on
              commandline) to find the latest directory with this path:
              <year>/dataset.week.#
              The program will find the latest <year> and within
              that year the dataset.week.# with highest #.  The output
              directories created will be put within this directory.
  
              Setting --customweekdir will cause program to use 'celppdir'
              path.
  
              Setting the --createweekdir flag will instruct this
              program to create a new directory for the current
              celpp week/year before invoking running any stage
              processing.
  
              NOTE: CELPP weeks start on Friday and end on Thursday
                    and week # follows ISO8601 rules so week numbers
                    at the end and start of the year are a bit
                    wonky.
  
              Breakdown of behavior of program is defined by
              value passed with --stage flag:
  
              If --stage 'import'
  
              In this stage 4 files are downloaded from urls specified
              by --compinchi and --pdbfileurl flags on the commandline
              into stage.1.dataimport directory.
  
              The tsv files are (--pdbfileurl flag sets url to
              download these files from):
  
              new_release_structure_nonpolymer.tsv
              new_release_structure_sequence.tsv
              new_release_crystallization_pH.tsv
  
              The ich file is (--compinchi flat sets url to
              download this file from):
  
              Components-inchi.ich
  
              If --stage 'blast'
  
              Verifies stage.1.dataimport exists and has 'complete'
              file.  Also the --blastdir path must exist and within a
              'current' symlink/directory must exist and within that a
              'complete' file must also reside. If both conditions
              are met then the 'blast' stage is run and output stored
              in stage.2.blastnfilter.  Requires --pdbdb to be set
              to a directory with valid PDB database files.

              If --stage 'proteinligprep'

              Verifies stage.2.blastnfilter exists and has 'complete'
              file.  If complete, this stage runs which invokes program
              set in --proteinligprep flag to prepare pdb and inchi files
              storing output in stage.3.proteinligprep.  --pdbdb flag
              must also be set when calling this stage.

              If --stage 'glide'

              Verifies stage3.proteinligprep exists and has a 'complete'
              file within it.  If complete, this stage runs which invokes
              program set in --glide flag to perform docking via glide
              storing output in stage.4.glide

              If --stage 'evaluation'

              Finds all stage.4.<algo> directories with 'complete' files
              in them which do not end in name 'webdata' and runs
              script set via --evaluation parameter storing the result of
              the script into stage.5.<algo>.evaluation. --pdbdb flag
              must also be set when calling this stage.
              
  
    positional arguments:
      celppdir              Base celpp directory
      
    optional arguments:
      -h, --help            show this help message and exit
      --blastdir BLASTDIR   Parent directory of blastdb. There should exist a
                            "current" symlink or directory that contains the db.
                            NOTE: Required parameter for blast stage
      --email EMAIL         Comma delimited list of email addresses
      --createweekdir       Create new celpp week directory before running stages
      --customweekdir       Use directory set in celppdir instead of looking for
                            latest weekdir. NOTE: --createweekdir will create a
                            dataset.week.# dir under celppdir
      --stage STAGE         Comma delimited list of stages to run. Valid STAGES =
                            {import, blast, proteinligprep, glide, evaluation}
      --blastnfilter BLASTNFILTER
                            Path to BlastnFilter script
      --postanalysis POSTANALYSIS
                            Path to PostAnalysis script
      --proteinligprep PROTEINLIGPREP
                            Path to proteinligprep script
      --glide GLIDE         Path to glide docking script
      --evaluation EVALUATION
                            Path to evaluation script
      --pdbdb PDBDB         Path to PDB database files
      --compinchi COMPINCHI
                            URL to download Components-inchi.ich file fortask
                            stage.1.compinchi
      --pdbfileurl PDBFILEURL
                            URL to download new_release_structure_nonpolymer.tsv
                            ,new_release_structure_sequence.tsv, and
                            new_release_crystallization_pH.tsv files for task
                            stage.1.dataimport
      --log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Set the logging level
      --smtp SMTP           Sets smtpserver to use
      --smtpport SMTPPORT   Sets smtp server port
      --version             show program's version number and exit

