#!/bin/sh

# download submissions from box

for u in NNNNNN_autodockvina NNNNNN_rdock-smallbox \
    NNNNNN_oefred-smallbox NNNNNN NNNNNN NNNNNN NNNNNN
do
    mkdir $u
    cd $u

    y=2017
    for week in $(seq 1 52)
    do
        echo $y $week $u
        curl -s -S -O --netrc-file ~/.curl ftps://ftp.box.com/celppweekly/usersubmissions/$u/celpp_week${week}_${y}_dockedresults_${u}.tar.gz
    done

    y=2018
    for week in $(seq 1 28)
    do
        echo $y $week $u
        curl -s -S -O --netrc-file ~/.curl ftps://ftp.box.com/celppweekly/usersubmissions/$u/celpp_week${week}_${y}_dockedresults_${u}.tar.gz
    done

    cd ..
done
