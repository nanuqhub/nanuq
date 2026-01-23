#!/bin/bash

NANUQ_EXE="../../../cfgs/generic/BLD/bin/nanuq.exe"

#XIOS_EXE="/home/laurent/src/xios2-trunk_oa3/bin/xios_server.exe"
XIOS_EXE="/opt/xios-trunk/bin/xios_server.exe"

#DIR_NC_IN="/data/laurent/INPUT_NANUQ_DISTRIB"
DIR_NC_IN="/SUMMER/DATA_MEOM/MEOM-OPENDAP/NANUQ/INPUT_NANUQ_DISTRIB"

#FATM_DIR="/SUMMER/SASIP/model-forcings/atmo_forcing/ERA5_Arctic" ; # path to ERA5 Arctic atmospheric forcing (on original ERA5 grid) mind: shape = 1440 x 264 !!!
FATM_DIR="/SUMMER/SASIP/model-forcings/atmo_forcing/ERA5_Arctic"

dir_in="${DIR_NC_IN}/EGL12"
dir_in_stdl="${dir_in}/standalone"
dir_in_bdy="${dir_in}/BDY"


YEAR=1997

for dr in "${dir_in}" "${dir_in_stdl}" "${dir_in_bdy}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

rm so*_bdy*.nc vo*_bdy*.nc ; # 3D ocean BDYs not needed here!

for dr in "${FATM_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*${YEAR}*.nc .
    ln -sf ${dr}/*$((YEAR-1))*.nc .
    ln -sf ${dr}/*$((YEAR+1))*.nc .
done

for exe in "${NANUQ_EXE}" "${XIOS_EXE}"; do
    ln -sf ${exe} .
done

