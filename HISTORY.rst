.. :changelog:

History
-------

0.1.0 (2015-06-30)
---------------------

* First release on PyPI.

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

