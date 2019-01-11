#!/bin/bash
#==============================================================================
#
# Script to build a modified version of the MAGMA library.
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
  local IS_CRAY_XK7 # OLCF Titan or Chester
  [[ -n "${CRAYOS_VERSION:-}" ]] && IS_CRAY_XK7="YES" || IS_CRAY_XK7="NO"
  local IS_IBM_AC922 # OLCF Summit or Peak
  [[ -n "${LSF_BINDIR:-}" ]] && IS_IBM_AC922="YES" || IS_IBM_AC922="NO"

  # WARNING: these module loads MUST match those in scripts/cmake.sh
  if [ $IS_CRAY_XK7 = YES ] ; then
    if [ "$PE_ENV" = "PGI" ] ; then
      module unload PrgEnv-pgi
    fi
    module load PrgEnv-gnu
    module load cudatoolkit
    module load acml
    cp ../make.inc.titan make.inc
  elif [ $IS_IBM_AC922 = YES ] ; then
    module load gcc/6.4.0
    module load cuda
    cp ../make.inc.summit make.inc
  else
    echo "Unknown platform." 1>&2
    exit 1
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
