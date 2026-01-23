#!/bin/bash
#
######################################################
#SBATCH -J <CASE>
#SBATCH --nodes=1
#SBATCH --ntasks=<NCORES2BOOK>
#SBATCH --account=<ACCOUNT>
#SBATCH --mem=<MEM>
#SBATCH --time=<TJOB>
#SBATCH --output out_<CONFCASE>_<CTI>_<cproc_shape>_%J.out
#SBATCH --error  err_<CONFCASE>_<CTI>_<cproc_shape>_%J.err
######################################################
#
INTEL_ONEAPI="/opt/intel/oneapi"
NCDF_INTEL="/workdir/cryodyn/chekkim/checkelmerice/Intel/install_elmerice/libelmerice/netcdf-4.7.2-install"
IC_MPI_V="2021.6.0"
IC_CMP_V="2022.1.0"
#
. ${INTEL_ONEAPI}/compiler/${IC_CMP_V}/env/vars.sh
. ${INTEL_ONEAPI}/mpi/${IC_MPI_V}/env/vars.sh
export LD_LIBRARY_PATH=${NCDF_INTEL}/lib:${LD_LIBRARY_PATH}
#
