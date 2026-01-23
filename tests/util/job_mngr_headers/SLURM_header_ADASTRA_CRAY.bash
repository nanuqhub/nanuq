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
module load PrgEnv-cray/8.5.0
module load gcc/13.2.0
#
# /opt/cray/libfabric/1.20.1/lib64   => 
#
H5NC_DIR="/lus/home/CT1/hmg2840/brodeau/opt/hdf5_netcdf4_cray_par"
export LD_LIBRARY_PATH=${H5NC_DIR}/lib:${LD_LIBRARY_PATH}

export OMP_NUM_THREADS=1
#
