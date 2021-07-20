#!/bin/bash

# Upload any remaining tar files to Box and remove locally.

if [ "$#" -ne 1 ]; then
    echo 'Must supply one directory.'
    exit 1
fi

week="$1"
cd $week
for f in $(ls *.tar); do
    echo $week/$f
    curl -s -S --netrc-file ~/.curl -T $f \
    	-ftp-create-dirs ftps://ftp.box.com/celppweekly/pdb_archive/$week/$f && \
        rm $f
done
cd ..
