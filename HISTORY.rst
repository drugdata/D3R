.. :changelog:

History
-------

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

