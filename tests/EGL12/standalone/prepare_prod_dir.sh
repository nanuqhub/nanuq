#!/bin/bash

NANUQ_EXE="../../../cfgs/generic/BLD/bin/nanuq.exe"

fsrc="../../paths_nanuq_data.bash" ; # path to file containing info relative to current host !
if [ -f ${fsrc} ]; then
    . ${fsrc}
else
    echo "I cannot find file: ${fsrc} !  :("
    exit
fi

dir_in="${DIR_NC_IN}/EGL12"
dir_in_stdl="${dir_in}/standalone"
dir_in_bdy="${dir_in}/BDY"


YEAR=1997

for dr in "${dir_in}" "${dir_in_stdl}" "${dir_in_bdy}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

rm so*_bdy*.nc vo*_bdy*.nc ; # 3D ocean BDYs not needed here!

if [ "${FATM_ERA5_DIR}" = "" ]; then echo " PROBLEM: variable FATM_ERA5_DIR must be set!"; exit; fi

for dr in "${FATM_ERA5_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*${YEAR}*.nc .
    ln -sf ${dr}/*$((YEAR-1))*.nc .
    ln -sf ${dr}/*$((YEAR+1))*.nc .
done

for exe in "${NANUQ_EXE}" "${XIOS_EXE}"; do
    ln -sf ${exe} .
done

