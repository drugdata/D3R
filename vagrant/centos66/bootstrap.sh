#!/usr/bin/env bash

echo "Installing base packages"
yum install -y git
yum install -y epel-release
yum install -y python-pip
yum install -y python-argparse
yum install -y python-lockfile
yum install -y python-psutil
yum install -y python-biopython
yum install -y python-virtualenv
yum install -y python-tox
yum install -y pylint
yum install -y python-coverage

echo "pip installing some packages"
pip install xlsxwriter
pip install ftpretty
pip install wheel
pip install flake8

echo "adding rdkit yum repo"
# install rdkit from yum repo
pushd /etc/yum.repos.d
wget https://copr.fedorainfracloud.org/coprs/giallu/rdkit/repo/epel-6/giallu-rdkit-epel-6.repo
popd

echo "installing rdkit"
yum install -y python-rdkit

yum install -y openbabel

echo "Installing ncbi blast pre-reqs"
# install ncbi blast
yum install -y perl-Archive-Tar
yum install -y perl-List-MoreUtils
echo "Installing ncbi blast.  May take a few minutes..."
rpm -ivh ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-1.x86_64.rpm

# install vina
wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz
tar -zxf autodock_vina_1_1_2_linux_x86.tgz
mv autodock_vina_1_1_2_linux_x86/bin/* /usr/local/bin/.

echo "Installing mgl tools"
# install mgltools
pushd /usr/local
wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
tar -zxf mgltools_x86_64Linux2_1.5.6.tar.gz
pushd mgltools_x86_64Linux2_1.5.6
./install.sh -d /usr/local/mgltools
popd
popd

echo "Install UCSF chimera"
# wget https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=linux_x86_64/chimera-1.10.2-linux_x86_64.bin -O chimera.bin
# chomd u+x chimera.bin
# ./chimera.bin

