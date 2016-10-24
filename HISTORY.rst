.. :changelog:

History
-------

1.6.1 (2016-10-24)
-------------------

* Added createchallenge stage which is NOT a stage, but a fake stage
  that runs the following stages: makedb,import,blast,challengedata. Issue #92

* Moved logic to setup logging handlers to start of celpprunner.py to remove
  no handlers found error for d3r.celpp.util. Issue #91

* Fixed bug where participant_list.csv could not be parsed if file was 
  written with carriage return delimiters instead of newlines. Issue #93
  
* Cleaned up CELPPade by updating documentation and simplifying variable names.

1.6.0 (2016-10-13)
-------------------

* Evaluation task now emails results of evaluation to external 
  submitter. Issues #49,#81

* Adjusted files uploaded to ftp server in EvaluationTask to 
  reflect changes in output from genchallengedata.py script.
  Issues #79,#80

* Added WebDavFileTransfer class to enable upload & download
  of files via WebDa for celppade tools. Issue #76

* Added tsv files and Components-inchi.ich files to list of 
  files uploaded to ftp by DataImportTask. Issue #78

* Updated challenge data package readme.txt to include documentation
  for hiTanimoto. Issue #75

* Added pdb_seqres.txt.gz to list of files uploaded to ftp by
  MakeBlastDBTask. Issue #77.

* Download canonical tsv file in data import stage. Issue #84

* Added --rdkitpython flag to celpprunner.py and modified
  code to pass it to chimera_proteinligprep.py. Issue #88

* Fixed bug where celpprunner would fail if evaluation 
  stage is rerun with completed evaluation tasks. Issue #87

* Improved documentation in RMSD.txt. Issues #82,#83

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
  removed and -ns flag added to preserve stereo information.

1.3.0 (2016-06-29)
---------------------

* Fixed bug #45 where autodock vina task was being incorrectly
  fed proteinligprep as input. Code now feeds it chimeraprep.

* #28 Data import stage waits for TSV files to be updated before
  downloading.  

* #8 celpprunner will now kill blastnfilter if it runs beyond
  time set via --blastnfiltertimeout flag.

* #37 Added external docking submission task which downloads
  external docked results so they can be evaluated the same
  way as the internal docking programs.

* #44 Added utility function to call external processes. To
  reduce redundancy in the code base.


1.2.0 (2016-06-03)
---------------------

* proteinligprep.py and chimera_proteinligprep.py
  has been updated to work with genchallenge stage output.

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
  the default upload directory for that challenge data package.

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

* First release on PyPI.

