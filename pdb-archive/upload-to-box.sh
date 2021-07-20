#!/bin/bash

# upload an entire PDB directory to Box

if [ "$#" -ne 1 ]; then
    echo 'Must supply one directory.'
    exit 1
fi

week="$1"
cd $week
for d in $(ls); do
    echo $week/$d
    tar cf $d.tar $d && \
        curl -s -S --netrc-file ~/.curl -T $d.tar \
            --retry 5 --retry-delay 5 \
            --ftp-create-dirs ftps://ftp.box.com/celppweekly/pdb_archive/$week/$d.tar && \
        rm $d.tar
done
cd ..

# if we fail to upload any tar files, call clean-to-box.sh
# to retry uploading them.
./cleanup-to-box.sh $week
