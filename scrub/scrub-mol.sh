#!/bin/sh

# scrub mol files

sed -i -e '1s/.*/ligand/' -e '2s/.*/ CelppScrubbed/' -e '3s/.*//' $*
