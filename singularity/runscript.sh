#!/bin/bash

# added cause local libstdc++.so.6 was being loaded instead of one shipped with openeye library
export LD_LIBRARY_PATH="/usr/lib/python2.7/site-packages/openeye/libs/python2.7-ucs4-linux-x64-g++4.x:$LD_LIBRARY_PATH"


exec /usr/bin/celpprunner.py "$@"
