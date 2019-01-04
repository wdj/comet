#!/bin/bash
#==============================================================================
#
# Script to build a modified version of the Magma library.
#
# Usage: type the following from the respective Magma root directory, e.g.,
# from magma_minproduct:
#
# ../make_magma.sh
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function do_make
{
  if [ -n "${CRAYOS_VERSION:-}" ] ; then
    IS_CRAY_XK7="YES"
  else
    IS_CRAY_XK7="NO"
  fi

  if [ $IS_CRAY_XK7 = YES ] ; then
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
    #module load cuda/9.1.85
    module load cuda
    cp ../make.inc.summit make.inc
  fi

  module list 2>&1

  time make lib -j8
}

#==============================================================================

function main
{
  do_make "$@" 2>&1 | tee out_make.txt
}

#==============================================================================

main "$@"

#==============================================================================
