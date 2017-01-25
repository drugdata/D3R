#!/usr/bin/env bash


echo "Installing git"
yum install -y git

echo "Cloning puppet-d3r"
git clone https://github.com/drugdata/puppet-d3r.git d3r

echo "Running puppet apply"
/opt/puppetlabs/bin/puppet apply --modulepath /home/vagrant -e "include d3r"
