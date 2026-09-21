#!/bin/bash

NANUQ_EXE="../../../../cfgs/generic/BLD/bin/nanuq.exe"

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

# Import reference config files:
for ff in "namelist_dom_ref" "namelist_ice_ref"; do
    ln -sf ../../../control_ref/${ff} .
done
for ff in "axis_def_nanuq.xml" "context_nanuq.xml" "domain_def_nanuq.xml" "field_def_nanuq.xml" "grid_def_nanuq.xml"; do
    ln -sf ../../../control_ref/xios/${ff} .
done

for exe in "${NANUQ_EXE}" "${XIOS_EXE}"; do
    ln -sf ${exe} .
done

