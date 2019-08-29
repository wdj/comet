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
  # Location of this script.
  local SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  . $SCRIPT_DIR/_platform_init.sh

  if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
    true # skip
  elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    cp ../make.inc.titan make.inc
    export GPU_TARGET=sm35
    export NV_SM=" -gencode arch=compute_35,code=sm_35"
    export NV_COMP="-gencode arch=compute_35,code=compute_35" MIN_ARCH=350
    MAKE_ARGS=""
  elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
    export CUDA_DIR="${CUDA_DIR:-$OLCF_CUDA_ROOT}"
    cp ../make.inc.summit make.inc
    export GPU_TARGET=sm70
    export NV_SM=" -gencode arch=compute_70,code=sm_70"
    export NV_COMP="-gencode arch=compute_70,code=compute_70" MIN_ARCH=350
    MAKE_ARGS=""
  elif [ $COMET_PLATFORM = DGX2 ] ; then
    cp ../make.inc.summit make.inc
    export CUDA_DIR=$HOME/cuda
    export GPU_TARGET=sm70
    export NV_SM=" -gencode arch=compute_70,code=sm_70"
    export NV_COMP="-gencode arch=compute_70,code=compute_70" MIN_ARCH=350
    MAKE_ARGS="CC=$HOME/.linuxbrew/bin/gcc-6 CXX=$HOME/.linuxbrew/bin/g++-6"
  elif [ $COMET_PLATFORM = GPUSYS2 ] ; then
    cp ../make.inc.summit make.inc
    export CUDA_DIR=/usr/local/cuda-10.1
    export GPU_TARGET=sm70
    export NV_SM=" -gencode arch=compute_75,code=sm_75"
    export NV_COMP="-gencode arch=compute_75,code=compute_75" MIN_ARCH=350
    MAKE_ARGS="CC=$(spack location --install-dir gcc)/bin/g++"
  else
    echo "Unknown platform." 1>&2
    exit 1
  fi

  time make lib $MAKE_ARGS -j8
}

#==============================================================================

function main
{
  do_make "$@" 2>&1 | tee out_make.txt
}

#==============================================================================

main "$@"

#==============================================================================
