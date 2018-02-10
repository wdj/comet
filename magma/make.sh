#!/bin/bash
#==============================================================================
#
# Script to build Magma library.
#
# Usage: type the following from the Magma root directory, e.g.,
# from magma_minproduct-1.6.2:
#
# ../make.sh
#
#==============================================================================

function do_make
{
  if [ -n "$CRAYOS_VERSION" ] ; then
    if [ "$PE_ENV" = "PGI" ] ; then
      module unload PrgEnv-pgi
    fi
    module load PrgEnv-gnu
    module load cudatoolkit
    module load acml
    cp ../make.inc.titan make.inc
  else #---IBM
    #if [ -n "$OLCF_XL_ROOT" ] ; then
    module load gcc
    #i
    module load cuda
    cp ../make.inc.summit make.inc
  fi

  time make lib -j8
}

#==============================================================================

do_make "$@" 2>&1 | tee out_make.txt

#==============================================================================
