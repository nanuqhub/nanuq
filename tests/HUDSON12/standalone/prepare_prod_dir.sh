#!/bin/bash

#DIR_NC_IN="/data/laurent/INPUT_NANUQ_DISTRIB"
#FATM_DIR="/data/laurent/HUDSON12/FATM" ; # path to atmospheric forcing interpolated on HUDSON12 domain

#DIR_NC_IN="/data/gcm_setup/INPUT_NANUQ_DISTRIB"
DIR_NC_IN="/home/laurent/tmp/INPUT_NANUQ_DISTRIB"

#FATM_DIR="/data/gcm_setup/HUDSON12/FATM_1997_HUDSON12" ; # path to atmospheric forcing interpolated on HUDSON12 domain
FATM_DIR="/home/laurent/tmp/FATM_1997_HUDSON12"


#if [ "${1}" = "" ]; then
#    echo "USAGE: ${0} <res_in_km>"
#    exit
#fi

dir_in="${DIR_NC_IN}/HUDSON12"
dir_in_stdl="${dir_in}/standalone"


for dr in "${dir_in}" "${dir_in_stdl}" "${FATM_DIR}"; do
    if [ ! -d ${dr} ];            then echo " PROBLEM: ${dr} does not exist!"; exit; fi
    ln -sf ${dr}/*.nc .
done

ln -sf ../../../cfgs/generic/BLD/bin/nanuq.exe .


