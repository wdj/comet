#!/bin/bash -l
#==============================================================================
#
# Script to build genomics metrics code.
#
#==============================================================================

if [ -n "$CRAYOS_VERSION" ] ; then
  # For Titan or Chester
  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  module load cudatoolkit
  module load acml
else #---IBM
  # For Summit or Peak
  module load spectrum-mpi/10.2.0.0-20180508 #FIX
  module load gcc/6.4.0
  #module load cuda/9.1.85
  module load cuda
fi

time make -j4 VERBOSE=1

if [ $? = 0 ] ; then
  time make install
  exit $?
else
  exit $?
fi

#==============================================================================
