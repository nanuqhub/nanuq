#!/bin/bash
################################
#SBATCH --account=hmg2840
#SBATCH --constraint=HPDA
#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --job-name=<CASE>
#SBATCH -o out_<CONFCASE>_<CTI>_<cproc_shape>_%J.out
#SBATCH -e err_<CONFCASE>_<CTI>_<cproc_shape>_%J.err
#SBATCH --time=<TJOB>
################################
