#!/bin/bash
################################
#SBATCH --account=<ACCOUNT>
#SBATCH --job-name=<CASE>
#SBATCH --constraint=MI250
#
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#
#SBATCH --mem=<MEM>
#
#SBATCH -o out_<CONFCASE>_<CTI>_<cproc_shape>_%J.out
#SBATCH -e err_<CONFCASE>_<CTI>_<cproc_shape>_%J.err
#SBATCH --time=<TJOB>
################################
#
module purge
export ROCM_PATH="/opt/rocm-7.2.0"
export OMPI_ROCM_PATH="/lus/home/CT1/hmg2840/brodeau/opt/openMPI5_rocm7p2"
export NCDF_ROCM_PATH="/lus/home/CT1/hmg2840/brodeau/opt/hdf5_netcdf4_rocm7p2_par"
export LD_LIBRARY_PATH=${ROCM_PATH}/lib:${OMPI_ROCM_PATH}/lib:${NCDF_ROCM_PATH}/lib:${LD_LIBRARY_PATH}
export PATH=${ROCM_PATH}/bin:${OMPI_ROCM_PATH}/bin:${NCDF_ROCM_PATH}/bin:${PATH}
#
export OMP_NUM_THREADS=1
#
