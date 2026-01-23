#!/bin/bash

#DIR_NC_IN="/DATA/IO/INPUT_NANUQ_DISTRIB"
DIR_NC_IN="/home/laurent/tmp/INPUT_NANUQ_DISTRIB"


if [ "${1}" = "" ]; then
    echo "USAGE: ${0} <res_in_km>"
    exit
fi

RESKM=${1}

lok="2 4 10"
if [ "`echo ${lok} | grep ${RESKM}`" = "" ]; then
    echo "Available resolutions in km are: ${lok} !"
    exit
fi
   

dir_in="${DIR_NC_IN}/CHANNELD17/${RESKM}km"

if [ ! -d ${dir_in} ]; then echo " PROBLEM: ${dir_in} does not exist!"; exit; fi

ln -sf ${dir_in}/*.nc .

ln -sf ../../cfgs/generic/BLD/bin/nanuq.exe .

# Picking the right namelists:
ln -sf namelist_dom_cfg.${RESKM}km namelist_dom_cfg
ln -sf namelist_ice_cfg.${RESKM}km namelist_ice_cfg

