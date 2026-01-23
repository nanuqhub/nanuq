#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
#
# ==========
# agrifpp.sh
# ==========
#
# ----------------------------
# Preform AGrif pre-processing
# ----------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ agrifpp.sh
#
#
# DESCRIPTION
# ===========
#
#
# Preprocess file using the conv in NEMOFILES directory
# Standard preprocessed files are stored in NEMOFILES/ppsrc/nemo
# Source files are stored under NEMOFILES/obj
# Include filess  in NEMOFILES/inc
# Note that agrif2model.F90 should not be preprocess (standard one) 
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./agrifpp.sh FILE_TO_PROCESS
# 
# TODO
# ====
#
# option debug
#
#
# EVOLUTIONS
# ==========
#
# $Id: agrifpp.sh 2143 2010-10-04 12:49:55Z rblod $
#
#
#
#   * creation
#
#-
MYDIR=$1
MYFILE=$(basename "$2")

if [ "$MYFILE" == "agrif2model.f90" ];then
   # generic case
   if [ -d ${MYDIR}/WORK ]; then
      \cp ${MYDIR}/WORK/${MYFILE/.f90/.F90} ${MYDIR}/NEMOFILES/obj/$MYFILE
   # DOMAINcfg case
   elif [ -d ${MYDIR}/src ]; then
      \cp ${MYDIR}/src/${MYFILE/.f90/.F90} ${MYDIR}/NEMOFILES/obj/$MYFILE
   fi
else
   cd ${MYDIR}/NEMOFILES/ppsrc/nemo
   ${MYDIR}/NEMOFILES/conv ${MYDIR}/NEMOFILES/agrif_oce.in -rm -incdir ${MYDIR}/NEMOFILES/inc -comdirout ${MYDIR}/NEMOFILES/obj -convfile ${MYFILE} > /dev/null
fi
