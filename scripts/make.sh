#!/bin/bash -l
#==============================================================================
#
# Script to build genomics metrics code.
#
#==============================================================================

if [ -n "$CRAYOS_VERSION" ] ; then
  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  module load cudatoolkit
  module load acml
else #---IBM
  #if [ -n "$OLCF_XL_ROOT" ] ; then
  module load gcc/6.4.0
  #fi
  module load cuda
fi

time make VERBOSE=1

if [ $? = 0 ] ; then
  time make install
  exit $?
else
  exit $?
fi

#==============================================================================
