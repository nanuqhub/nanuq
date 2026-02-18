#!/bin/bash

NEMO_EXE="/home/laurent/NEMO/NEMOv4.2.2/cfgs/HUDSON12_OPA_OA3/BLD/bin/nemo.exe"

NANUQ_EXE="../../../cfgs/generic_cpl_oce/BLD/bin/nanuq.exe"

fsrc="../../paths_nanuq_data.bash" ; # path to file containing info relative to current host !
if [ -f ${fsrc} ]; then
    . ${fsrc}
else
    echo "I cannot find file: ${fsrc} !  :("
    exit
fi

dir_in="${DIR_NC_IN}/HUDSON12"
dir_in_cpl="${dir_in}/cpl_OCE"

for dr in "${dir_in}" "${dir_in_cpl}" "${ERA5_HUDSON12_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

for exe in "${NANUQ_EXE}" "${NEMO_EXE}"; do
    ln -sf ${exe} .
done
