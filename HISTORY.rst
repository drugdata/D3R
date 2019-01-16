.. :changelog:

History
-------

1.11.2 (2019-01-16)
--------------------
* Fix for failed test in 1.11.1

1.11.1 (2019-01-16)
--------------------

* Email participant when submission is empty.
  `Issue #197 <https://github.com/drugdata/D3R/issues/197>`_

* Allow spaces in passwords in ftp config file.
  `Issue #198 <https://github.com/drugdata/D3R/issues/198>`_


1.11.0 (2018-09-07)
--------------------

* Fixed minor bug where unit tests for molfilevalidator.py failed when valid
  openeye license was available.
  `Issue #177 <https://github.com/drugdata/D3R/issues/177>`_

* Added retry to external submission download to prevent complete failure of
  download in event of network hiccup
  `Issue #181 <https://github.com/drugdata/D3R/issues/181>`_

* Fixed bug where median value incorrectly calculated in evaluate.py. Instead
  of averaging middle two values for case of even number of elements, the old
  code just chose the latter value.
  `Issue #183 <https://github.com/drugdata/D3R/issues/183>`_

* Fixed bug where median value incorrectly calculated in post_evaluation.py.
  Instead of averaging middle two values, the old code just chose the latter
  value.
  `Issue #184 <https://github.com/drugdata/D3R/issues/184>`_

* Updated BlastNFilterSummary class (used by celppreports.py) to calculate
  number of targets found by counting target .txt files if value is not in
  summary.txt file
  `Issue #187 <https://github.com/drugdata/D3R/issues/187>`_

* Added call to website REST service in EvaluationTask and PostEvaluationTask
  to persist evaluation results to website. This change added a new flag
  --websiteserviceconfig that requires a configuration file.
  `Issue #188 <https://github.com/drugdata/D3R/issues/188>`_

* Fixed bug in evaluate.py wait_and_check() function which caused it to 
  prematurely give up on an alignment.
  `Issue #189 <https://github.com/drugdata/D3R/issues/189>`_

* Changed evaluate.py so it now writes final.log directly to output directory.
  `Issue #190 <https://github.com/drugdata/D3R/issues/190>`_

* Moved code under if __main__ into main() function to facilitate testing.
  `Issue #191 <https://github.com/drugdata/D3R/issues/191>`_

* Add RMSD.csv to files uploaded in EvaluationTask get_uploadable_files().
  `Issue #192 <https://github.com/drugdata/D3R/issues/192>`_

* Changed logging in evaluate.py to use module logger instead of root logger
  `Issue #194 <https://github.com/drugdata/D3R/issues/194>`_

1.10.0 (2018-01-09) 
--------------------

* **Breaking change** Added new flag to celpprunner.py --smtpconfig
  that lets a user specify smtp configuration information. 
  At the same time removed --replyto --smtp, and --smtpport since
  those options would be confusing to have also in place.
  `Issue #166 <https://github.com/drugdata/D3R/issues/166>`_ 

* molfilevalidator.py has new command line option --excludedir to
  let caller ignore directories when looking for mol files.
  `Issue #167 <https://github.com/drugdata/D3R/issues/167>`_

* Fixed bug where a network hiccup raised an exception during
  DataImportTask causing celpprunner.py to exit prematurely.
  `Issue #168 <https://github.com/drugdata/D3R/issues/168>`_

* Improved error message output when a participant uploads a
  malformed challenge data package. 
  `Issue #169 <https://github.com/drugdata/D3R/issues/169>`_

* Version of D3R is now output in post evaluation summary email.
  `Issue #172 <https://github.com/drugdata/D3R/issues/172>`_

* Fixed bug where celpprunner.py was not sending an email to
  people in --summaryemail list if a task fails.
  `Issue #171 <https://github.com/drugdata/D3R/issues/171>`_

* Post evaluation summary email now outputs lines with NA
  values for any submissions that failed
  `Issue #130 <https://github.com/drugdata/D3R/issues/130>`_

1.9.2 (2017-10-30)
--------------------

* Added molfilevalidator.py to validate D3R submission tarfiles.
  `Issue #165 <https://github.com/drugdata/D3R/issues/165>`_

* evaluation.py now generates RMSD.json which is a JSON version of
  RMSD files. `Issue #143 <https://github.com/drugdata/D3R/issues/143>`_

* Added a fix to deal with out of memory errors encountered in 
  blastnfilter stage. `Issue #5 <https://github.com/drugdata/D3R/issues/5>`_

1.9.1 (2017-08-21)
--------------------

* Fixed bug where Apo targets not getting pocket center correctly defined
  `Issue #151 <https://github.com/drugdata/D3R/issues/151>`_ 

1.9.0 (2017-06-23)
--------------------

* EvaluationTask now records evaluate.py task exit code in a file
  `Issue #134 <https://github.com/drugdata/D3R/issues/134>`_

* Symmetry filter added to blastnfilter.py
  `Issue #145 <https://github.com/drugdata/D3R/issues/145>`_

* In RMSD.txt,RMSD.csv renamed Medium to Median and swapped values
  for Maximum and Minimum
  `Issue #144 <https://github.com/drugdata/D3R/issues/144>`_

* Added note about new values in parenthesis in individual results
  email
  `Issue #142 <https://github.com/drugdata/D3R/issues/142>`_

1.8.0 (2017-05-18)
--------------------

* EvaluationTask modified to pass path to blastnfilter task to evaluate.py
  `Issue #139 <https://github.com/drugdata/D3R/issues/139>`_

* Median RMSD added to post_evaluation.py outputs
  `Issue #136 <https://github.com/drugdata/D3R/issues/136>`_

* Added histogram of RMSD scores to post_evaluation.py outputs
  `Issue #137 <https://github.com/drugdata/D3R/issues/137>`_

* Continuously output the analysis result into the pickle csv and txt files
  `Issue #133 <https://github.com/drugdata/D3R/issues/133>`_

* Align the docked complex using the binding site alignment for each of the crystal template and calculate the RMSD, if the binding site alignment failed, then the whole protein alignment will be applied

* Improve the extraction step in the evaluate.py to ensure the RMSD calculating was only applied to the docked ligand but not all others ligand like solvents or co-factors

* Add the ligand center calculation step in the evaluate.py to output the distance of the docked ligand with the crystal ligand. Also calculate the distance between the original LMCSS ligand center with the crystal ligand center

* Update the genchallengedata.py to extract Apo chain which is closed to the LMCSS ligand `Issue #135 <https://github.com/drugdata/D3R/issues/135>`_


1.7.3 (2017-04-11)
--------------------

* requirements.txt and setup.py modified to require biopython 
  v1.6.8, which is the last version compatible with python2.6


1.7.2 (2017-04-11)
--------------------

* Bug fix, evaluation of a submission will not fail if one
  candidate in that submission fails.


1.7.1 (2017-03-27)
--------------------

* Bug fix, evaluation task will not fail if call to external
  script fails. `Issue #129 <https://github.com/drugdata/D3R/issues/129>`_

1.7.0 (2017-03-21)
--------------------

* Added Post Evaluation Stage which summarizes evaluations
  of all docking submissions sending email to addresses set
  by --summaryemail flag
  `Issue #111 <https://github.com/drugdata/D3R/issues/111>`_
  `Issue #110 <https://github.com/drugdata/D3R/issues/110>`_
  `Issue #100 <https://github.com/drugdata/D3R/issues/100>`_

* If any task/stage fails report that via email to addresses
  set by --summaryemail flag
  `Issue #125 <https://github.com/drugdata/D3R/issues/125>`_ 

* Order of evaluation of external dock submissions can now be 
  dictated by a new column in participant_list.csv
  `Issue #124 <https://github.com/drugdata/D3R/issues/124>`_

* A timeout has now been added for blastnfilter stage to prevent
  the task/stage from running too long. Default is 24 hours, but
  can be changed with --blastnfiltertimeout flag.
  `Issue #8 <https://github.com/drugdata/D3R/issues/8>`_

* A timeout has now been added for evaluation stage to prevent
  the task/stage from running too long. Default is 24 hours, but
  can be changed with --evaluationtimeout flag.
  `Issue #123 <https://github.com/drugdata/D3R/issues/123>`_ 

* Fixed bug in blastnfilter to correctly rank the list of 
  hiTanimoto candidates.

1.6.8 (2017-03-07)
------------------

* Celpprunner puts lock file within specific week directory. 
  `Issue #122 <https://github.com/drugdata/D3R/issues/122>`_

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

