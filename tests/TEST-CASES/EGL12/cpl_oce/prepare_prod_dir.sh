#!/bin/bash

NEMO_EXE="/home/laurent/NEMO/NEMOv4.2.2/cfgs/EGL12_OPA_OA3/BLD/bin/nemo.exe"

NANUQ_EXE="../../../../cfgs/generic_cpl_oce/BLD/bin/nanuq.exe"

fsrc="../../paths_nanuq_data.bash" ; # path to file containing info relative to current host !
if [ -f ${fsrc} ]; then
    . ${fsrc}
else
    echo "I cannot find file: ${fsrc} !  :("
    exit
fi

dir_in="${DIR_NC_IN}/EGL12"
dir_in_cpl="${dir_in}/cpl_OCE"
dir_in_bdy="${dir_in}/BDY"


YEAR=1997

for dr in "${dir_in}" "${dir_in_cpl}" "${dir_in_bdy}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

if [ "${FATM_ERA5_DIR}" = "" ]; then echo " PROBLEM: variable FATM_ERA5_DIR must be set!"; exit; fi

for dr in "${FATM_ERA5_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*${YEAR}*.nc .
    ln -sf ${dr}/*$((YEAR-1))*.nc .
    ln -sf ${dr}/*$((YEAR+1))*.nc .
done

# Import reference config files:
for ff in "namelist_dom_ref" "namelist_ice_ref"; do
    ln -sf ../../../control_ref/${ff} .
done
for ff in "axis_def_nanuq.xml" "context_nanuq.xml" "domain_def_nanuq.xml" "field_def_nanuq.xml" "grid_def_nanuq.xml"; do
    ln -sf ../../../control_ref/xios/${ff} .
done

for exe in "${NANUQ_EXE}" "${NEMO_EXE}"  "${XIOS_OA3_EXE}"; do
    if [ ! -f ${exe} ]; then echo "ERROR: \`${exe}\` does not exist!"; exit; fi
    ln -sf ${exe} .
done
