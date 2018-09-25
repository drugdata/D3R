#!/bin/bash -l

module load python/1
module load singularity

export SCHRODINGER=/opt/schrodinger
export SCHROD_LICENSE_FILE=@@SHROD_LICENSE_FILE@@


celppdir="/oasis/tscc/scratch/@@USER@@/celpp/data"
cd $celppdir


singularity run --bind /oasis --bind /proc --bind /state @@D3R_IMAGE@@ --stage extsubmission  --log DEBUG --email @@EMAIL_LIST@@ --ftpconfig @@BOX_CONFIG@@ --smtpconfig @@SMTP_CONFIG@@ $celppdir
ecode=$?

echo "Done time is `date` and exit code is: $ecode"

exit $ecode

