===============================
Drug Design Data Resource CELPP Runner
===============================

.. image:: https://img.shields.io/travis/nbcrrolls/D3R.svg
        :target: https://travis-ci.org/nbcrrolls/D3R

.. image:: https://img.shields.io/pypi/v/D3R.svg
        :target: https://pypi.python.org/pypi/D3R


Drug Design Data Resource is a suite of software to enable 
filtering, docking, and scoring of new sequences from wwpdb.

* Free software: BSD license
* Documentation: https://D3R.readthedocs.org.

Features
--------

 * Works with Python 2.6, 2.7, 3.0, 3.2, 3.4

Installation
------------

Pip install is coming, but in the meantime:

.. code:: bash

  git clone https://github.com/nbcrrolls/D3R
  cd D3R
  make install

Usage
-----

Run

.. code:: bash
  
  celpprunner.py -h
  usage: celpprunner.py [-h] [--blastdir BLASTDIR] [--email EMAIL] --stage
                      {blast,dock,score} --blastnfilter BLASTNFILTER
                      [--log {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                      [--smtp SMTP] [--smtpport SMTPPORT] [--version]
                      celppdir

              Runs last 3 stages (blast, dock, & score) of CELPP
              processing pipeline (http://www.drugdesigndata.org)

              CELPP processing pipeline relies on a set of directories
              with specific structure. The pipeline runs a set of stages
              Each stage has a numerical value and a name. The numerical
              value denotes order and the stage name identifies separate
              tasks to run in the stage.

              The filesystem structure of the stage is:

              stage.<stage number>.<task name>

              Only 1 stage is run per invocation of this program and the
              stage to be run is defined via the required --stage flag.

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

              Regardless of the stage specified, this program will
              examine the 'celppdir' (last argument passed on
              commandline) to find the latest directory with this path:
              <year>/dataset.week.#
              The program will find the latest <year> and within
              that year the dataset.week.# with highest #.  The output
              directories created will be put within this directory.

              Breakdown of behavior of program is defined by
              value passed with --stage flag:

              If --stage 'blast'

              Verifies stage.1.dataimport exists and has 'complete'
              file.  Also the --blastdir path must exist and within a
              'current' symlink/directory must exist and within that a
              'complete' file must also reside. If both conditions
              are met then the 'blast' stage is run and output stored
              in stage.2.blastnfilter

              If --stage 'dock'

              Verifies stage2.blastnfilter exists and has a 'complete'
              file within it.  If complete, this program will run fred
              docking and store output in stage.3.fred.  As new
              algorithms are incorporated additional stage.3.<algo> will
              be created and run.

              If --stage 'score'

              Finds all stage.3.<algo> directories with 'complete' files
              in them and invokes appropriate scoring algorithm storing
              results in stage.4.<algo>.scoring.
              

  positional arguments:
    celppdir              Base celpp directory

  optional arguments:
    -h, --help            show this help message and exit
    --blastdir BLASTDIR   Parent directory of blastdb. There should exist a
                          "current" symlink or directory that contains the db.
    --email EMAIL         Comma delimited list of email addresses
    --stage {blast,dock,score}
                          Stage to run blast = blastnfilter (2), dock = fred &
                          other docking algorithms (3), score = scoring (4)
    --blastnfilter BLASTNFILTER
                          Path to BlastnFilter script
    --log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                          Set the logging level
    --smtp SMTP           Sets smtpserver to use
    --smtpport SMTPPORT   Sets smtp server port
    --version             show program's version number and exit


* TODO
