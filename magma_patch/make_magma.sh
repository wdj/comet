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
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  cp ../$COMET_MAGMA_MAKE_INC make.inc
  local _CMGA_=$COMET_MAGMA_GPU_ARCH
  env CUDA_DIR=$CUDA_ROOT \
    GPU_TARGET=sm$_CMGA_ MIN_ARCH=350 \
    NV_SM=" -gencode arch=compute_${_CMGA_},code=sm_${_CMGA_}" \
    NV_COMP=" -gencode arch=compute_${_CMGA_},code=compute_${_CMGA_}" \
  time make lib CC=$COMET_CXX_SERIAL_COMPILER -j8
}

#==============================================================================

function main
{
  do_make "$@" 2>&1 | tee out_make.txt
}

#==============================================================================

main "$@"

#==============================================================================
