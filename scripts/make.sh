#!/bin/bash -l
#==============================================================================
#
# Script to build genomics metrics code.
#
#==============================================================================

if [ "$PE_ENV" = "PGI" ] ; then
  module swap PrgEnv-pgi PrgEnv-gnu
fi
module load cudatoolkit
module load acml

make VERBOSE=1

if [ $? = 0 ] ; then
  make install
fi

#==============================================================================
