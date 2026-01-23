#!/bin/bash

#DIR_NC_IN="/DATA/IO/INPUT_NANUQ_DISTRIB"
DIR_NC_IN="/home/laurent/tmp/INPUT_NANUQ_DISTRIB"


if [ "${1}" = "" ]; then
    echo "USAGE: ${0} <res_in_km>"
    exit
fi

RESKM=${1}

lok="2 10"
if [ "`echo ${lok} | grep ${RESKM}`" = "" ]; then
    echo "Available resolutions in km are: ${lok} !"
    exit
fi

case ${RESKM} in
    "2")  DT="240"      ; # NANUQ time step in seconds
            NTS="32400" ; # number of time steps to go...
            ;;
    "10") DT="1200"     ; # NANUQ time step in seconds
            NTS="6480"  ; # number of time steps to go...
            ;;
esac

dir_in="${DIR_NC_IN}/ROTATION/${RESKM}km"

if [ ! -d ${dir_in} ]; then echo " PROBLEM: ${dir_in} does not exist!"; exit; fi

ln -sf ${dir_in}/*.nc .

ln -sf ../../cfgs/generic/BLD/bin/nanuq.exe .

# Building the right namelists:
sed -e s/"<RESKM>"/"${RESKM}km"/g -e s/"<DT>"/"${DT}"/g -e s/"<NTS>"/"${NTS}"/g namelist_dom_cfg.tmplt > namelist_dom_cfg
sed -e s/"<RESKM>"/"${RESKM}km"/g                                               namelist_ice_cfg.tmplt > namelist_ice_cfg
