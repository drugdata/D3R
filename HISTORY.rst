.. :changelog:

History
-------

1.6.7 (2017-02-07)
------------------

* Increased retry time for challengedata download to 1 hour

* Updates to evaluation suffix structure (Allows for more complex suffixes after guids)

* Fixed bugs in CELPPade error messages

1.6.6 (2017-01-23)
------------------

* Add retry to external submission downloads. 
  `Issue #112 <https://github.com/drugdata/D3R/issues/112>`_

1.6.5 (2017-01-13)
------------------

* Fix the bug to skip invalid holo hit instead of terminating the whole blastnfilter process. 
  `Issue #89 <https://github.com/drugdata/D3R/issues/89>`_

* Fix the bug in the set sequence fuction where it original complains about the lower cased chain ID

1.6.4 (2017-01-10)
------------------

* Remove intermediate files with pymol prefix from challenge data
  package. `Issue #73 <https://github.com/drugdata/D3R/issues/73>`_

* Place docked files in top-level target directory in submissions.
  `Issue #86 <https://github.com/drugdata/D3R/issues/86>`_

* Fixed sporadically failing unit test. `Issue #104 <https://github.com/drugdata/D3R/issues/104>`_

* Add evaluation chain permuter implementation and tests. `Issue #107 <https://github.com/drugdata/D3R/issues/107>`_

1.6.3 (2016-12-01)
-------------------

* Modified ParticipantDatabase to optionally strip off _# from guid
  when doing search for Participant. This is to handle case where
  single participant has multiple submissions to CELPP. `Issue #98 <https://github.com/drugdata/D3R/issues/98>`_

* Switched os.getlogin() calls to  pwd.getpwuid(os.getuid())[0] 
  cause os.getlogin() was raising OSError on Travis. `Issue #102 <https://github.com/drugdata/D3R/issues/102>`_

* Fixed bug where evaluation task completed email had log messages
  from other evaluation task emails in them. `Issue #99 <https://github.com/drugdata/D3R/issues/99>`_

* Updated readme.txt file in challenge data package to reflect
  use of new_release_structure_sequence_canonical.tsv instead of
  new_release_structure_sequence.tsv file. `Issue #97 <https://github.com/drugdata/D3R/issues/97>`_

1.6.2 (2016-10-26)
-------------------

* Fixed bug where large amounts of output to standard out/err caused
  celpprunner.py to exit due to an exception from smtplib due to 
  very large email. `Issue #95 <https://github.com/drugdata/D3R/issues/95>`_

1.6.1 (2016-10-24)
-------------------

* Added createchallenge stage which is NOT a stage, but a fake stage
  that runs the following stages: makedb,import,blast,challengedata. `Issue #92 <https://github.com/drugdata/D3R/issues/92>`_

* Moved logic to setup logging handlers to start of celpprunner.py to remove
  no handlers found error for d3r.celpp.util. `Issue #91 <https://github.com/drugdata/D3R/issues/91>`_

* Fixed bug where participant_list.csv could not be parsed if file was 
  written with carriage return delimiters instead of newlines. `Issue #93 <https://github.com/drugdata/D3R/issues/93>`_
  
* Cleaned up CELPPade by updating documentation and simplifying variable names

* Version of d3r is now written to 'start' file in each stage/task. `Issue #94 <https://github.com/drugdata/D3R/issues/94>`_

1.6.0 (2016-10-13)
-------------------

* Evaluation task now emails results of evaluation to external 
  submitter. Issues `#49 <https://github.com/drugdata/D3R/issues/49>`_ , `#81 <https://github.com/drugdata/D3R/issues/81>`_

* Adjusted files uploaded to ftp server in EvaluationTask to 
  reflect changes in output from genchallengedata.py script.
  Issues `#79 <https://github.com/drugdata/D3R/issues/79>`_ , `#80 <https://github.com/drugdata/D3R/issues/80>`_

* Added WebDavFileTransfer class to enable upload & download
  of files via WebDa for celppade tools. `Issue #76 <https://github.com/drugdata/D3R/issues/76>`_ 

* Added tsv files and Components-inchi.ich files to list of 
  files uploaded to ftp by DataImportTask. `Issue #78 <https://github.com/drugdata/D3R/issues/78>`_

* Updated challenge data package readme.txt to include documentation
  for hiTanimoto. `Issue #75 <https://github.com/drugdata/D3R/issues/75>`_

* Added pdb_seqres.txt.gz to list of files uploaded to ftp by
  MakeBlastDBTask. `Issue #77 <https://github.com/drugdata/D3R/issues/77>`_

* Download canonical tsv file in data import stage. `Issue #84 <https://github.com/drugdata/D3R/issues/84>`_

* Added --rdkitpython flag to celpprunner.py and modified
  code to pass it to chimera_proteinligprep.py. `Issue #88 <https://github.com/drugdata/D3R/issues/88>`_

* Fixed bug where celpprunner would fail if evaluation 
  stage is rerun with completed evaluation tasks. `Issue #87 <https://github.com/drugdata/D3R/issues/87>`_

* Improved documentation in RMSD.txt. Issues `#82 <https://github.com/drugdata/D3R/issues/82>`_ , `#83 <https://github.com/drugdata/D3R/issues/83>`_

1.5.0 (2016-09-11)
--------------------

* Modified blastnfilter candidate txt file by adding hiTanimoto and 
  adding more information to hiResHolo and SMCSS.

* Genchallengedata.py modified to keep single chains for all holo
  proteins (LMCSS, SMCSS, hiResHolo, hiTanimoto)

* In blastnfilter, hiResHolo now only has top structure 
  reported and only one chain. Where top structure is 
  highest resolution hit.

1.4.0 (2016-08-11)
--------------------

* Fixed issue #66 Change candidate category names. Largest is now LMCSS,
  Smallest is now SMCSS, Apo is now HiResApo, Holo is now HiResHolo

1.3.4
--------------------

* Fixed issue #58 in chimera_proteinligprep.py code now uses rdkit 
  for 3d conf gen instead of babel

1.3.3 (2016-07-18)
--------------------

* Fixed issue #60 where challenge data package was NOT being
  uploaded to remote server

1.3.2 (2016-07-12)
--------------------

* Removed #8 blastnfilter timeout since it was causing blastnfilter
  script to hang.

* Blastnfilter.py now uses argparse to parse command line arguments

* Added loggging support into blastnfilter.py 

1.3.1 (2016-07-01)
---------------------

* In proteinligprep.py ligprep command modified. -s 1 -g flags 
  removed and -ns flag added to preserve stereo information

1.3.0 (2016-06-29)
---------------------

* Fixed bug #45 where autodock vina task was being incorrectly
  fed proteinligprep as input. Code now feeds it chimeraprep

* #28 Data import stage waits for TSV files to be updated before
  downloading

* #8 celpprunner will now kill blastnfilter if it runs beyond
  time set via --blastnfiltertimeout flag

* #37 Added external docking submission task which downloads
  external docked results so they can be evaluated the same
  way as the internal docking programs

* #44 Added utility function to call external processes. To
  reduce redundancy in the code base


1.2.0 (2016-06-03)
---------------------

* proteinligprep.py and chimera_proteinligprep.py
  has been updated to work with genchallenge stage output

* vinadocking.py and glidedocking.py now output receptor as pdb
  and ligand as mol.

* evaluate.py modified to accept new output format as described
  here:  https://github.com/drugdata/D3R/wiki/Proposed-challenge-docked-results-file-structure

* ProteinLigPrepTask #41 modified to use ChallengeDataTask as input

* ChimeraProteinLigPrepTask #30 modified to use ChallengeDataTask as input

1.1.0 (2016-05-24)
---------------------

* ChallengeDataTask now uploads challenge data package 
  (celpp_week##_##.tar.gz) to 'challengedata' directory on
  ftp if ftpconfig is set properly.  This is in addition, to
  the default upload directory for that challenge data package

* Added a header line in readme.txt of challenge data package
  to denote start of Blastnfilter summary output.

* Not part of production release, but added prototype vagrant 
  configuration to enable easy creation of a VM that can run
  celpprunner.

1.0.0 (2016-05-12)
---------------------

* Added chimeraprep stage to prepare data with Chimera 
  (issue #32)
 
* Added challengedata stage to generate challenge data package (issue #22)
  and added genchallengedata.py script which does the work (issue #21)

* Added vina stage to run docking with autodock vina (issue #15)
  and added vinadocking.py script to run the docking

* Modified D3rTask to write error message into 'error' file (issue #12)

* Added celppreports.py to provide summary reports (issue #14)

* Modified DataImportTask to compare entries in tsv file with 
  data in pdb_seqres.txt in makeblastdb stage.  As part of this
  fix made dataimport stage dependent on makeblastdb stage so
  the order is now stage.1.makeblastdb => stage.2.dataimport =>
  stage.3.blastnfilter... (issue #16)

0.1.0 (2015-06-30)
---------------------

* First release on PyPI

