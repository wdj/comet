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

function script_dir
{
  local TARGET_FILE="$0"

  if [ $(uname -s) != "Darwin" ] ; then
    echo $(dirname "$(readlink -f "$TARGET_FILE")")
    return
  fi

  cd $(dirname $TARGET_FILE)
  TARGET_FILE=$(basename $TARGET_FILE)
  # Iterate down a (possible) chain of symlinks
  while [ -L "$TARGET_FILE" ] ; do
    TARGET_FILE=$(readlink $TARGET_FILE)
    cd $(dirname $TARGET_FILE)
    TARGET_FILE=$(basename $TARGET_FILE)
  done
  # Compute the canonicalized name by finding the physical path
  # for the directory we're in and appending the target file.
  local PHYS_DIR=$(pwd -P)
  local RESULT=$PHYS_DIR/$TARGET_FILE
  echo $(dirname $RESULT)
}

#==============================================================================

function do_make
{
  # Location of this script.
  local SCRIPT_DIR=$(script_dir)
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  if [ ${USE_CUDA:-OFF} = ON ] ; then
    cp ../$COMET_MAGMA_MAKE_INC make.inc
    local _CMGA_=$COMET_MAGMA_GPU_ARCH
    time env CUDA_DIR=$CUDA_ROOT \
      GPU_TARGET=sm$_CMGA_ MIN_ARCH=350 \
      NV_SM=" -gencode arch=compute_${_CMGA_},code=sm_${_CMGA_}" \
      NV_COMP=" -gencode arch=compute_${_CMGA_},code=compute_${_CMGA_}" \
    make lib CC=$COMET_CXX_SERIAL_COMPILER -j8
  fi

  if [ ${USE_HIP:-OFF} = ON ] ; then
    cp make.inc-examples/make.inc.hip_openblas make.inc

    sed -i -e 's/GPU_TARGET = gfx803 gfx900 gfx901/GPU_TARGET = gfx906 gfx908/' make.inc
    sed -i -e 's/DEVCCFLAGS  = -O3 -DNDEBUG -DADD_/DEVCCFLAGS  = -O3 -DNDEBUG -DADD_ --amdgpu-target=gfx906,gfx908/' make.inc
#    if [ ${USE_BLIS:-OFF} = ON ] ; then
#      sed -i -e 's/lopenblas/lblis/' make.inc
#      export OPENBLASDIR=$PWD/../../blis/blis
    if [ ${USE_LAPACK:-OFF} = ON ] ; then
      sed -i -e 's/lopenblas/lrefblas/' make.inc
      sed -i -e 's/-frecursive/-frecursive -fPIC/' make.inc
      export OPENBLASDIR=$PWD/../../lapack/lapack
    else
      sed -i -e 's/lopenblas/lsci_cray/' make.inc
      export OPENBLASDIR=$CRAY_LIBSCI_PREFIX
    fi

    env HIPDIR=$HIP_PATH make -f make.gen.hipMAGMA_*

    # tools/codegen.py fails on non-ascii characters in files, thus:
    sed -i -e '122d' magma_*blas_hip/zlarfg.hip.cpp
    sed -i -e '93d' magma_*blas_hip/zlarfg-v2.hip.cpp
    sed -i -e '100d' magma_*blas_hip/zlarfgx-v2.hip.cpp
    sed -i -e '128d' magma_*blas_hip/zlarfgx-v2.hip.cpp
    sed -i -e '595d' testing/magma_*_generate.cpp
    sed -i -e '13d' sparse/blas/zgeellrtmv.cu
    sed -i -e '115d' sparse/blas/zgeellrtmv.cu
    sed -i -e '64d' sparse/blas/zgeellrtmv.cu
    sed -i -e 's/.*#include .*cublas.h.*//' include/*.h

    env HIPDIR=$HIP_PATH  make lib -j16
  fi
}

#==============================================================================

function main
{
  do_make "$@" 2>&1 | tee out_make.txt
}

#==============================================================================

main "$@"

#==============================================================================
