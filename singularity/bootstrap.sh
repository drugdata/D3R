#!/usr/bin/env bash

echo "Installing base packages"
yum install -y git epel-release 
yum install -y wget tar unzip zip gcc xauth squashfs-tools libarchive-devel
# yum install -y python-pip 
# yum install -y python-argparse python-lockfile python-psutil
#yum install -y python-biopython python-virtualenv python-tox
#yum install -y pylint python-coverage libXft mesa-* openbabel
#yum install -y perl-Archive-Tar perl-List-MoreUtils xauth pymol

#VERSION=2.5.0
#wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
#tar xvf singularity-$VERSION.tar.gz
#cd singularity-$VERSION
#./configure --prefix=/usr/local
#make
#make install

wget https://github.com/singularityware/singularity/releases/download/2.3.2/singularity-2.3.2.tar.gz
tar xvf singularity-2.3.2.tar.gz
cd singularity-2.3.2
./configure --prefix=/usr/local
make
make install

