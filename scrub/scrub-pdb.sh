#!/bin/sh

# scrub pdb files

sed -i '/^\(ATOM\|CONECT\|END\|ENDMDL\|HETATM\|MASTER\)/!d' $*
