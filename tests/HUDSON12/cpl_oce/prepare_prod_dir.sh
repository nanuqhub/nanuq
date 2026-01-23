#!/bin/bash

NEMO_EXE="/home/laurent/NEMO/NEMOv4.2.2_BBM/cfgs/HUDSON12_OPA_OA3/BLD/bin/nemo.exe"
NANUQ_EXE="../../../cfgs/generic_cpl_oce/BLD/bin/nanuq.exe"

DIR_NC_IN="/data/laurent/INPUT_NANUQ_DISTRIB"
FATM_DIR="/data/laurent/HUDSON12/FATM_1997_HUDSON12" ; # path to atmospheric forcing interpolated on HUDSON12 domain

dir_in="${DIR_NC_IN}/HUDSON12"
dir_in_cpl="${dir_in}/cpl_OCE"

for dr in "${dir_in}" "${dir_in_cpl}" "${FATM_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

for exe in "${NANUQ_EXE}" "${NEMO_EXE}"; do
    ln -sf ${exe} .
done
