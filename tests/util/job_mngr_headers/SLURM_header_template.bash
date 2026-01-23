#!/bin/bash
################################
#SBATCH -N <NNODES_TOT>
#SBATCH -n <NCORES2BOOK>
#SBATCH -J <CASE>
#SBATCH -o out_<CONFCASE>_<CTI>_<cproc_shape>_%J.out
#SBATCH -e err_<CONFCASE>_<CTI>_<cproc_shape>_%J.err
#SBATCH --time=<TJOB>
#SBATCH <QUEUE_OPTION>
