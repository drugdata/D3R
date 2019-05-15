#!/bin/sh

# copy RMSD.csv from results into appropriate finals/ directory.


process() {
    local u=$1
    local week=$2
    local y=$3

    fdir="final/${u}/celpp_week${week}_${y}_dockedresults_${u}"
    if [ -d $fdir ]
    then
        rfile="results/${y}_week${week}/${u}/RMSD.csv"

        if [ -f $rfile ]
        then
            #echo mv $rfile ${fdir}/
            mv $rfile ${fdir}/
        #else
            #echo not $rfile
        fi
    fi
}


for u in NNNNN_autodockvina NNNNN_rdock-smallbox \
    NNNNN_oefred-smallbox NNNNNN NNNNNN NNNNNN NNNNNN \
    NNNNN_glide
do
    y=2017
    for week in $(seq 1 52)
    do
        process $u $week $y
    done

    y=2018
    for week in $(seq 1 28)
    do
        process $u $week $y
    done

done
