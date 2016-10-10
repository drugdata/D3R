#!/bin/bash


# download pdb
echo "Be sure to run pdb_updater.sh before going to next step"
echo "Hit Ctrl-C to cancel... Sleeping 10 seconds"

sleep 10

#  run celpprunner
celpprunner.py --createweekdir --stage makedb,import,blast,challengedata,proteinligprep,chimeraprep,glide,vina  --pdbdb /vagrant/pdb.extracted --log DEBUG --compinchi http://ligand-expo.rcsb.org/dictionaries --proteinligprep proteinligprep.py --glide glidedocking.py --skipimportwait /vagrant/celpp
