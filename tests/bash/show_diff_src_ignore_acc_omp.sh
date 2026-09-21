#!/usr/bin/env bash


src_dir_1=${1}
src_dir_2=${2}


diff -r -B -I '!$omp' -I '!$acc' -I '!Nothing for' -I '_OPENACC' -I '_OPENMP'  ${src_dir_1}   ${src_dir_2}






