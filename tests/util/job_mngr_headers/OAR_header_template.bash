#!/bin/bash
################################
#OAR -n <CASE>
#OAR -l /nodes=<NNODES_TOT>/core=<NCORES_TOT>,walltime=<TJOB>
#OAR -O out_<CONFCASE>_<CTI>_<cproc_shape>_%jobid%.out
#OAR -E err_<CONFCASE>_<CTI>_<cproc_shape>_%jobid%.err
#OAR --project <QUEUE>
###OAR -p "cpumodel = '<cpu_type>'"
################################
#
ulimit -s unlimited
#
source ${HOME}/.nix-profile/bin/iccvars.sh -arch intel64 -platform linux
#

