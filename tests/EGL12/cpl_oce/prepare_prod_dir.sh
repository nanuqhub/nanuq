#!/bin/bash

NEMO_EXE="/home/laurent/NEMO/NEMOv4.2.2/cfgs/EGL12_OPA_OA3/BLD/bin/nemo.exe"

NANUQ_EXE="../../../cfgs/generic_cpl_oce/BLD/bin/nanuq.exe"

fsrc="../../paths_nanuq_data.bash" ; # path to file containing info relative to current host !
if [ -f ${fsrc} ]; then
    . ${fsrc}
else
    echo "I cannot find file: ${fsrc} !  :("
    exit
fi

DIR_NC_IN="/data/laurent/INPUT_NANUQ_DISTRIB"
FATM_DIR="/SUMMER/SASIP/model-forcings/atmo_forcing/ERA5_Arctic" ; # path to ERA5 Arctic atmospheric forcing (on original ERA5 grid) mind: shape = 1440 x 264 !!!

dir_in="${DIR_NC_IN}/EGL12"
dir_in_cpl="${dir_in}/cpl_OCE"
dir_in_bdy="${dir_in}/BDY"


YEAR=1997

for dr in "${dir_in}" "${dir_in_cpl}" "${dir_in_bdy}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

for dr in "${FATM_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*${YEAR}*.nc .
    ln -sf ${dr}/*$((YEAR-1))*.nc .
    ln -sf ${dr}/*$((YEAR+1))*.nc .
done


for exe in "${NANUQ_EXE}" "${NEMO_EXE}"  "${XIOS_OA3_EXE}"; do
    ln -sf ${exe} .
done
