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
    # Build for Titan or Chester
    if [ "$PE_ENV" = "PGI" ] ; then
      module unload PrgEnv-pgi
    fi
    module load PrgEnv-gnu
    module load cudatoolkit
    module load acml
    cp ../make.inc.titan make.inc
  else #---IBM
    # Build for Summit or Peak
    module load gcc/6.4.0
    module load cuda/9.1.85
    cp ../make.inc.summit make.inc
  fi

  time make lib -j8
}

#==============================================================================

do_make "$@" 2>&1 | tee out_make.txt

#==============================================================================
