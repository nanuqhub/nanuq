#!/bin/bash


host=`hostname | cut -d '.' -f2`
if [ "`hostname | cut -d '.' -f2-3`" = "u-ga.fr" ]; then host=`hostname | cut -d '.' -f1`; fi

# Defaults:
export NBCPN=0 ; # will trigger error
export NP_TC=4 ; # default number of procs in // to use for test-cases...
export NEMO_REPO_ROOT="${HOME}/NEMO"
export JOBMNGR="none"
export QUEUE_OPTION="none"
export QUEUE_EXTRA="none"

export CP="rsync -LavP" ; # Default command to copy a file
export LNK="ln -sf"     ; # Default command for symbolic linking

export I_RST_DIR_TMP=0 ; # 0 => RST dir to use, containing generic restarts is found in the normal `<CONF>-I/RST/` dir ...
#                         # 1 => it is found in the production aka TMPDIR

case ${host} in
    "luitel"|"ige-mcpc-36"|"mcp-oceannext-01") export ARCH="LUITEL"
                                               export NBCPN=4
                                               export DATA_DIR="/data/gcm_setup"
                                               export FATM_DIR="${DATA_DIR}/FATM"
                                               export DIR_STOR_READ_ROOT="${DATA_DIR}"
                                               export DIR_STOR_WRIT_ROOT="/scratch/work"
                                               export DIR_STOR_SAVE_ROOT="/scratch/work"                                               
                                               ;;
    "ige-mcpc-39")    export ARCH="POMME"
                      export DATA_DIR="/opt/data"
                      export FATM_DIR="${DATA_DIR}/FATM"
                      ;;
    "frazilo"|"ige-osugb1-s115") export ARCH="FRAZILO"
                                 export NBCPN=30
                                 export DATA_DIR="/data/laurent"
                                 export FATM_DIR="${DATA_DIR}/atmo_forcing"
                                 export DIR_STOR_READ_ROOT="${DATA_DIR}"
                                 export DIR_STOR_WRIT_ROOT="/home/data/laurent/tmp"
                                 export DIR_STOR_SAVE_ROOT="/home/data/laurent/tmp"
                                 export DIR_FATM_ROOT=""
                                 export XIOS_HOME="/opt/xios-trunk"
                                 export NP_TC=16 ; # default number of procs in // to use for test-cases...
                                 ;;
    "merlat"        ) export ARCH="MERLAT"
                      export NBCPN=8
                      export DATA_DIR="/DATA/IO"
                      export FATM_DIR="${DATA_DIR}/FATM"
                      export DIR_STOR_READ_ROOT="${DATA_DIR}"
                      export DIR_STOR_WRIT_ROOT="${DATA_DIR}/tmp"
                      export DIR_STOR_SAVE_ROOT="${DATA_DIR}/tmp"
                      export DIR_FATM_ROOT=""
                      export XIOS_HOME="/opt/xios-trunk"
                      ;;
    "ige-meom-cal1" ) export ARCH="MEOMCAL1"
                      export DATA_DIR="/mnt/meom/workdir/brodeau"
                      #export FATM_DIR="${DATA_DIR}/FATM"
                      export FATM_DIR="/mnt/summer/SASIP/model-forcings/atmo_forcing"
                      export NCKS="ncks -h -a"
                      #
                      export INTEL_ONEAPI="/mnt/meom/workdir/brodeau/opt/intel/oneapi"
                      export NCDF_INTEL="/mnt/meom/workdir/brodeau/opt/hdf5_netcdf4_intel_par"
                      #
                      export PATH=${INTEL_ONEAPI}/compiler/latest/linux/bin/intel64:${PATH}
                      export LD_LIBRARY_PATH=${INTEL_ONEAPI}/compiler/latest/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
                      export CPATH=${INTEL_ONEAPI}/compiler/include:${INTEL_ONEAPI}/compiler/include/intel64:${CPATH}
                      #
                      export INTEL_MPI_DIR="${INTEL_ONEAPI}/mpi/latest"
                      export PATH=${INTEL_MPI_DIR}/bin:${PATH}
                      export LD_LIBRARY_PATH=${INTEL_MPI_DIR}/lib:${LD_LIBRARY_PATH}
                      export CPATH=${INTEL_MPI_DIR}/include:${CPATH}
                      #
                      export LD_LIBRARY_PATH=${NCDF_INTEL}/lib:${LD_LIBRARY_PATH}
                      ;;
    "login1"|"login2"|"login3"|"login4"|"login5"|"login6"|"login7" ) export ARCH="ADASTRA"
                                                   export NBCPN=192 ; # number of cores per node =  MACHINE SPECIFIC !!!
                                                   #export NBCPN=384 ; # number of cores per node =  MACHINE SPECIFIC !!!
                                                   export JOBMNGR="SLURM"
                                                   #export QUEUE_OPTION="--account=c1512020"
                                                   export QUEUE_OPTION="--account=hmg2840"
                                                   export QUEUE_EXTRA="--constraint=GENOA"
                                                   export DATA_DIR="/lus/store/CT1/hmg2840/brodeau" ; # STOREDIR
                                                   export FATM_DIR="${DATA_DIR}/atmo_forcing"
                                                   export NEMO_REPO_ROOT="/lus/store/CT1/hmg2840/brodeau/NEMO"
                                                   export DIR_STOR_READ_ROOT="/lus/store/CT1/hmg2840/brodeau" ; # STOREDIR
                                                   #export DIR_STOR_WRIT_ROOT=${SCRATCHDIR}
                                                   #export DIR_STOR_SAVE_ROOT=${SCRATCHDIR}
                                                   export DIR_STOR_WRIT_ROOT=/lus/scratch/CT1/hmg2840/brodeau
                                                   export DIR_STOR_SAVE_ROOT=/lus/scratch/CT1/hmg2840/brodeau
                                                   export I_RST_DIR_TMP=1
                                                   export LNK=${CP}     ; # no links !!!
                                                   #export DIR_STOR_SAVE_ROOT=${WORKDIR}/tmp
                                                   ;;
    "lacroix"       ) export ARCH="LACROIX"
                      export DATA_DIR="/data2/work"
                      export FATM_DIR="${DATA_DIR}/FATM"
                      ;;
    "fram" )          export ARCH="FRAM"
                      export DATA_DIR="/cluster/projects/nn9878k/brodeau"
                      export FATM_DIR="${DATA_DIR}/FATM"
                      ;;
    "jackzilla" )     export ARCH="JACKZILLA"
                      export NBCPN=33
                      export DATA_DIR="/data1/nobackup/laurent"
                      export FATM_DIR="${DATA_DIR}/atmo_forcing"
                      export NCKS="ncks -h --no-alphabetize"
                      export DIR_STOR_READ_ROOT=/data1/nobackup/laurent
                      export DIR_STOR_WRIT_ROOT=/data1/nobackup/laurent/tmp
                      export DIR_STOR_SAVE_ROOT=/data2/nobackup/laurent/tmp
                      export NP_TC=16 ; # default number of procs in // to use for test-cases...
                      ;;
    "f-dahu"       )  export ARCH="DAHU"
                      export DATA_DIR="/bettik/brodeaul/GCMIO"
                      ;;
    *               ) echo "Unknow architecture: ${host}"
                      exit
                      ;;
esac

#
#export DIR_STOR_READ_ROOT="${DATA_DIR}"
#export DIR_STOR_WRIT_ROOT="${DATA_DIR}/tmp"
#export DIR_STOR_SAVE_ROOT="${DATA_DIR}/tmp"


echo
echo " *** Based on the HOST: ${host} => ARCH = ${ARCH}"; echo
echo "      - DATA_DIR = ${DATA_DIR}"
echo "      - FATM_DIR = ${FATM_DIR}"
echo ""
echo "      - DIR_STOR_READ_ROOT = ${DIR_STOR_READ_ROOT}"
echo "      - DIR_STOR_WRIT_ROOT = ${DIR_STOR_WRIT_ROOT}"
echo "      - DIR_STOR_SAVE_ROOT = ${DIR_STOR_SAVE_ROOT}"
echo
