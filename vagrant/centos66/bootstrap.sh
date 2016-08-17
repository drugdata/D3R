#!/usr/bin/env bash

echo "Installing base packages"
yum install -y git epel-release
yum install -y python-pip
yum install -y python-argparse python-lockfile python-psutil
yum install -y python-biopython python-virtualenv python-tox
yum install -y pylint python-coverage libXft mesa-* openbabel
yum install -y perl-Archive-Tar perl-List-MoreUtils xauth

echo "pip installing some packages"
pip install xlsxwriter
pip install ftpretty
pip install wheel
pip install flake8
pip install lockfile --upgrade

#echo "adding rdkit yum repo"
# install rdkit from yum repo
#pushd /etc/yum.repos.d
#wget https://copr.fedorainfracloud.org/coprs/giallu/rdkit/repo/epel-6/giallu-rdkit-epel-6.repo
#popd

#echo "installing rdkit"
#yum install -y python-rdkit

echo "Installing ncbi blast pre-reqs"
# install ncbi blast
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

#
# install rdkit
pushd /vagrant
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod a+x Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh -b -p /opt/miniconda2
export PATH="/opt/miniconda2/bin:$PATH"
conda update --yes conda
conda install -y -c rdkit rdkit=2016.03.3

if [ -f "/vagrant/chimera-1.10.2-linux_x86_64.bin" ] ; then
  pushd /vagrant
  chmod a+x chimera-1.10.2-linux_x86_64.bin -o /opt/chimera -y
  echo "/opt/chimera" | ./chimera-1.10.2-linux_x86_64.bin 
  popd
else 
  echo "MANUAL TASK -- Install UCSF chimera since it wasn't found in /vagrant"
  echo "1) Download chimera by visiting: https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=linux_x86_64/chimera-1.10.2-linux_x86_64.bin"
  echo "   In a browser, accept license and download binary."
  echo "   Copy binary to directory with VagrantFile & bootstrap.sh file"
  echo ""
  echo "2) Connect to vm via vagrant ssh"
  echo "3) Become root sudo -u root /bin/bash"
  echo "4) cd /vagrant"
  echo "5) Run chmod a+x chimera-1.10.2-linux_x86_64.bin"
  echo "6) Run ./chimera-1.10.2-linux_x86_64.bin"
  echo "   and follow instructions using defaults"
  echo ""
fi


if [ -f "/vagrant/Schrodinger_Suites_2016-2_Linux-x86_64.tar" ] ; then
  pushd /vagrant
  tar -xf Schrodinger_Suites_2016-2_Linux-x86_64.tar
  pushd /vagrant/Schrodinger_Suites_2016-2_Linux-x86_64
  ./INSTALL -d `pwd` -b -s /vagrant/schrodinger -k /usr/tmp -t /vagrant/schrodinger/thirdparty mmshare*.gz glide*.gz maestro*gz
  popd
  popd
else
  echo "MANUAL INSTALL -- /vagrant/Schrodinger_Suites_2016-2_Linux-x86_64.tar not found"
  echo ""
  echo "Download Schrodinger_Suites_2016-2_Linux-x86_64.tar"
  echo "tar -xf Schrodinger_Suites_2016-2_Linux-x86_64.tar"
  echo "cd Schrodinger_Suites_2016-2_Linux-x86_64"
  echo "./INSTALL"
fi
