#!/bin/bash
################################
#SBATCH --account=<ACCOUNT>
#SBATCH --job-name=<CASE>
#SBATCH --constraint=GENOA
#SBATCH --nodes=<NNODES_TOT>
#SBATCH --ntasks=<NCORES2BOOK>
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH -o out_<CONFCASE>_<CTI>_<cproc_shape>_%J.out
#SBATCH -e err_<CONFCASE>_<CTI>_<cproc_shape>_%J.err
#SBATCH --time=<TJOB>
################################
#
export INTEL_ONEAPI="/lus/home/softs/intel/oneapi"
export IC_MPI_V="2021.6.0"
export IC_CMP_V="2022.1.0"
export COMP_LIB_DIR=${INTEL_ONEAPI}/compiler/${IC_CMP_V}/linux/compiler/lib/intel64_lin
. ${INTEL_ONEAPI}/compiler/${IC_CMP_V}/env/vars.sh
. ${INTEL_ONEAPI}/mpi/${IC_MPI_V}/env/vars.sh
export H5NC_DIR="/lus/home/CT1/hmg2840/brodeau/opt/hdf5_netcdf4_intel_par"
#
export LD_LIBRARY_PATH="${H5NC_DIR}/lib:${LD_LIBRARY_PATH}"
#
unset I_MPI_SHM_SEND_TINY_MEMCPY_THRESHOLD
unset I_MPI_DAPL_DIRECT_COPY_THRESHOLD
#
export OMP_NUM_THREADS=1
#
