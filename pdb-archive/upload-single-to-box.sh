#!/bin/bash

# upload a single file to Box

if [ "$#" -ne 1 ]; then
    echo 'Must supply one file.'
    exit 1
fi

f="$1"

echo $f

curl -s -S --netrc-file ~/.curl -T $f \
    --retry 20 --retry-delay 10 \
    --ftp-create-dirs ftps://ftp.box.com/celppweekly/pdb_archive/$f
