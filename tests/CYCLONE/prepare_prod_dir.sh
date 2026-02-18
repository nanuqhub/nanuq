#!/bin/bash

if [ "${1}" = "" ]; then
    echo "USAGE: ${0} <res_in_km>"
    exit
fi

RESKM=${1}

fsrc="../paths_nanuq_data.bash" ; # path to file containing info relative to current host !
if [ -f ${fsrc} ]; then
    . ${fsrc}
else
    echo "I cannot find file: ${fsrc} !  :("
    exit
fi



lok="2 4 10"
if [ "`echo ${lok} | grep ${RESKM}`" = "" ]; then
    echo "Available resolutions in km are: ${lok} !"
    exit
fi

case ${RESKM} in
    "2")  DT="120"     ; # NANUQ time step in seconds
            NTS="1440"   ; # number of time steps to go...
            ;;
    "4")  DT="240"     ; # NANUQ time step in seconds
            NTS="720"   ; # number of time steps to go...            
            ;;
    "10") DT="600"    ; # NANUQ time step in seconds
            NTS="288"   ; # number of time steps to go...
            ;;
esac

dir_in="${DIR_NC_IN}/CYCLONE/${RESKM}km"

if [ ! -d ${dir_in} ]; then echo " PROBLEM: ${dir_in} does not exist!"; exit; fi

ln -sf ${dir_in}/*.nc .

ln -sf ../../cfgs/generic/BLD/bin/nanuq.exe .

# Building the right namelists:
sed -e s/"<RESKM>"/"${RESKM}km"/g -e s/"<DT>"/"${DT}"/g -e s/"<NTS>"/"${NTS}"/g namelist_dom_cfg.tmplt > namelist_dom_cfg
sed -e s/"<RESKM>"/"${RESKM}km"/g                                               namelist_ice_cfg.tmplt > namelist_ice_cfg
