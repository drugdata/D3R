#!/bin/sh

# anonymize user names

new=$(echo $* | sed -e 's/NNNNNN_autodockvina/autodockvina/' \
    -e 's/NNNNNN_rdock-smallbox/rdock/' \
    -e 's/NNNNNN_oefred-smallbox/oefred/' \
    -e 's/NNNNNN_glide/glide/' \
    -e 's/NNNNNN/external_participant1/' \
    -e 's/NNNNNN/external_participant2/' \
    -e 's/NNNNNN/external_participant4/' \
    -e 's/NNNNNN/external_participant3/' )

echo mv $* $new
#mv $* $new
