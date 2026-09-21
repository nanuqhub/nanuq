#!/usr/bin/env bash

#set +H

NP_PAR=24 ; # Number of jobs in parallel...

spath=`dirname ${0}`

#ioverwrite=0
#if [ "${1}" = "-o" ]; then
#    echo; echo " * Will overwrite all *.F90 files..."
#    sleep 3
#    ioverwrite=1
#fi


files_with_acc=$(grep -rl --include="*90" '!$acc end data' src/)

echo


ijob=0

for ff in ${files_with_acc}; do
    ijob=$(( ijob + 1 ))
    echo ; echo
    echo "###############################################################"
    echo ${ff}
    echo "###############################################################"


    cat ${ff} | grep '!$acc end data'
    #
    echo "###############################################################"
    echo

    if [ $((ijob % NP_PAR)) -eq 0 ]; then
        wait
    fi    
    #
done

