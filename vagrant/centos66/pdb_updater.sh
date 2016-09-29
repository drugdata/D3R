#!/bin/bash

rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/pdb/ /vagrant/pdb
if [ $? != 0 ]; then
        echo -e "RSYNC failed rsync.rcsb.org::ftp_data/structures/divided/pdb/\n\nSincerely,\n$0"
        exit 1
fi

rm -rf /vagrant/pdb.extracted

cp -r /vagrant/pdb /vagrant/pdb.new
cd /vagrant/pdb.new

find . -name *.gz -exec gunzip {} \;
if [ $? != 0 ]; then
        echo -e "GUNZIP failed rsync.rcsb.org::ftp_data/structures/divided/pdb/\n\nSincerely,\n$0"
        exit 1
fi

cd /vagrant
mv pdb.new pdb.extracted

