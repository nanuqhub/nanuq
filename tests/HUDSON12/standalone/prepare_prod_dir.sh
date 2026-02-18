#!/bin/bash

NANUQ_EXE="../../../cfgs/generic/BLD/bin/nanuq.exe"

fsrc="../../paths_nanuq_data.bash" ; # path to file containing info relative to current host !
if [ -f ${fsrc} ]; then
    . ${fsrc}
else
    echo "I cannot find file: ${fsrc} !  :("
    exit
fi

dir_in="${DIR_NC_IN}/HUDSON12"
dir_in_stdl="${dir_in}/standalone"


for dr in "${dir_in}" "${dir_in_stdl}" "${ERA5_HUDSON12_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

ln -sf ${NANUQ_EXE} .


