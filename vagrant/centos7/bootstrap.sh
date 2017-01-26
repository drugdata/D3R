#!/usr/bin/env bash


echo "Installing git and cloning puppet-d3r"
yum install -y git

git clone https://github.com/drugdata/puppet-d3r.git

echo "Running puppet apply"
/opt/puppetlabs/bin/puppet apply --modulepath /home/vagrant/ --catalog /home/vagrant/puppet-d3r/metadata.json