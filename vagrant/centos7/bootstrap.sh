#!/usr/bin/env bash


echo "Installing git and cloning puppet-d3r"
yum install -y git

git clone https://github.com/drugdata/puppet-d3r.git

# append execution command cause I cannot for the life of
# me get puppet to use the catalog 
# TODO: fix this hack at some point

echo -e "\n# Run the d3r open class\nclass { 'd3r::open': }\n" >> /home/vagrant/puppet-d3r/manifests/open.pp

echo "Running puppet apply"
/opt/puppetlabs/bin/puppet apply /home/vagrant/puppet-d3r/manifests/open.pp
