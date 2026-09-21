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



lok="2 10"
if [ "`echo ${lok} | grep ${RESKM}`" = "" ]; then
    echo "Available resolutions in km are: ${lok} !"
    exit
fi

NDAYS_RUN=60 ; # number of model days to complete


case ${RESKM} in
    "2")  DT="240"      ; # NANUQ time step in seconds
            ;;
    "10") DT="1200"     ; # NANUQ time step in seconds
            ;;
esac

NTS=$((3600*24/DT*NDAYS_RUN)) ; # number of time steps to go...

dir_in="${DIR_NC_IN}/ROTATION/${RESKM}km"

# Import reference config files:
for ff in "namelist_dom_ref" "namelist_ice_ref"; do
    ln -sf ../../control_ref/${ff} .
done
for ff in "axis_def_nanuq.xml" "context_nanuq.xml" "domain_def_nanuq.xml" "field_def_nanuq.xml" "grid_def_nanuq.xml"; do
    ln -sf ../../control_ref/xios/${ff} .
done

# Domain, forcings, etc:
if [ ! -d ${dir_in} ]; then echo " PROBLEM: ${dir_in} does not exist!"; exit; fi
ln -sf ${dir_in}/*.nc .

# The exec!
ln -sf ../../../cfgs/generic/BLD/bin/nanuq.exe .

# Building the right namelists:
sed -e s/"<RESKM>"/"${RESKM}km"/g -e s/"<DT>"/"${DT}"/g -e s/"<NTS>"/"${NTS}"/g namelist_dom_cfg.tmplt > namelist_dom_cfg
sed -e s/"<RESKM>"/"${RESKM}km"/g                                               namelist_ice_cfg.tmplt > namelist_ice_cfg
