#!/usr/bin/env bash

#set +H

#NP_PAR=24 ; # Number of jobs in parallel...

CONV_SCRIPT="OACCtoOMP.sh"                  ; # name of script that converts OpenACC -> OpenMP, one file at a time
#                                          ; # => must be in the same directory as current script...
IN_SRC_DIR="./src"                         ; # directory containing sources to treat
MY_SRC_DIR="./cfgs/generic/MY_SRC"         ; # directory in which converted files will be saved

spath=`dirname ${0}`

if [ "${2}" = "" ]; then
    echo; echo " USAGE: $0 <PATH_DIR_INPUT_SRC> <PATH_DIR_TARGET_SRC>"
    exit
fi

IN_SRC_DIR="${1}"
MY_SRC_DIR="${2}"


if [ "${NP_PAR}" = "" ] || [ ${NP_PAR} -lt 1 ]; then
    echo " * Please set the variable \`NP_PAR\` to something >= 1 !"
    exit
fi

if [ ! -f ${spath}/${CONV_SCRIPT} ]; then
    echo " * PROBLEM: could not find scrip ${spath}/${CONV_SCRIPT} "; exit
fi

mkdir -p ${MY_SRC_DIR}

if [ "${VERSION_OMP}" != "4.x" ] && [ "${VERSION_OMP}" != "5.1" ] && [ "${VERSION_OMP}" != "cray" ]; then
    echo; echo "Please set the environment variable VERSION_OMP (4.x or 5.1 or cray) !"; echo
    echo "  Using a compiler that supports OMP version 5.1 or higher is highly recommended!"
    echo
    exit
fi

# Populating files that contain OpenACC directives:
files_with_acc=$(grep -rl --include="*90" '!$acc' ${IN_SRC_DIR}/)

echo
echo " The following files will be converted:"
for ff in ${files_with_acc}; do
    echo " * ${ff}"
done
echo
echo "  ==> will convert by batches of ${NP_PAR} files in parallel (variable NP_PAR) !"
echo
sleep 3
echo

ijob=0

for ff in ${files_with_acc}; do
    #
    ijob=$(( ijob + 1 ))
    #
    echo ; echo
    echo "###############################################################"
    echo ${ff}
    echo "###############################################################"
    
    fo=$( basename ${ff} ) ; # name of OMP-compliant file to be generated
    
    # Check if we have a version of this file specially adapted to OMP (mostly due to ROCm bugs...)
    if [ -f ${ff}.omp ]; then
        echo
        echo " * I found a this file: ${ff}.omp !"
        echo "   ==> and I will convert this one rather than ${ff}... "
        ff=${ff}.omp
        echo
    fi

    # 1 Back-up:
    #rsync -avP ${ff} ${ff}.openACC

    
    echo
    echo "   ==> will convert ${ff} to ${MY_SRC_DIR}/${fo} ! "
    #
    ${spath}/${CONV_SCRIPT} ${ff} ${MY_SRC_DIR} ${fo} > conversion_`basename ${ff}`.log &
    #
    echo "###############################################################"
    echo

    if [ $((ijob % NP_PAR)) -eq 0 ]; then
        echo "... waiting ..."
        wait
    fi
    #
done

wait

#if [ ${ioverwrite} -eq 1 ]; then
#    for ff in ${files_with_acc}; do
#        mv -f ${ff}.openMP ${ff}
#    done
#fi


echo
