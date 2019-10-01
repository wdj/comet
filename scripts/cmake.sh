#!/bin/bash -l
#==============================================================================
#
# Configure script for cmake.
#
# Usage: to configure and build code, perform the following steps:
#
# cd genomics_gpu/..
# mkdir build
# cd build
# ../genomics_gpu/scripts/cmake.sh # configure
# ../genomics_gpu/scripts/make.sh # make
#
# Relevant input variables:
#
# INSTALL_DIR - installation directory (default: see below)
# USE_MPI - ON or OFF (default), build with true MPI or for single process use
# BUILD_TYPE - Debug (default) or Release
# FP_PRECISION - SINGLE or DOUBLE (default) precision for floating point
#    operations (does not affect use of tensor cores for CCC)
# TESTING - ON or OFF (default), to build testing code
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function get_modules_used_magma
{
  if [ -n "$MODULEPATH" ] ; then
    (module -t list) 2>&1 | sort | egrep 'gcc|cuda' 
  fi
}

#==============================================================================

function main
{
  # Location of this script.
  local SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  #----------------------------------------------------------------------------
  #---Basic platform-specific settings.

  local CC_serial=g++

  if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
    true # skip
  #--------------------
  elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    local cc=$(which cc)
    local CC=$(which CC)
  #--------------------
  elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
    local cc=$(which mpicc)
    local CC=$(which mpiCC)
  #--------------------
  elif [ $COMET_PLATFORM = DGX2 ] ; then
    local cc=$HOME/.linuxbrew/bin/gcc-6
    local CC=$HOME/.linuxbrew/bin/g++-6
  #--------------------
  elif [ $COMET_PLATFORM = GPUSYS2 ] ; then
    local cc=$(spack location --install-dir gcc)/bin/gcc
    local CC=$(spack location --install-dir gcc)/bin/g++
    local CC_serial=$CC
  #--------------------
  elif [ $COMET_PLATFORM = EDISON ] ; then
    local cc=$(which cc)
    local CC=$(which CC)
  #--------------------
  elif [ $COMET_PLATFORM = LYRA ] ; then
    local cc=$(which gcc)
    local CC=$(which g++)
    CC_serial=hipcc
  #--------------------
  elif [ $COMET_PLATFORM = AMDINTERNAL ] ; then
    local cc=$(which gcc)
    local CC=$(which g++)
    CC_serial=hipcc
  #--------------------
  else
    echo "${0##*/}: Unknown platform." 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Cleanup of old files.

  rm -rf CMakeCache.txt CMakeFiles

  #----------------------------------------------------------------------------
  #---Set build type.

  [[ -z "${BUILD_TYPE:-}" ]] && local BUILD_TYPE=Debug #=Release
  if [ "$BUILD_TYPE" != "Debug" -a "$BUILD_TYPE" != "Release" ] ; then
    echo "${0##*/}: Invalid setting for BUILD_TYPE. $BUILD_TYPE" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  # Set whether to build with MPI stub library.

  [[ -z "${USE_MPI:-}" ]] && local USE_MPI=ON #=OFF
  if [ "$USE_MPI" != "ON" -a "$USE_MPI" != "OFF" ] ; then
    echo "${0##*/}: Invalid setting for USE_MPI. $USE_MPI" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Set floating point precision for certain calculations.

  [[ -z "${FP_PRECISION:-}" ]] && local FP_PRECISION=DOUBLE #=SINGLE
  if [ "$FP_PRECISION" != "SINGLE" -a "$FP_PRECISION" != "DOUBLE" ] ; then
    echo "${0##*/}: Invalid setting for FP_PRECISION. $FP_PRECISION" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Set whether to build unit test code.

  [[ -z "${TESTING:-}" ]] && local TESTING=OFF #=ON
  if [ "$TESTING" != "ON" -a "$TESTING" != "OFF" ] ; then
    echo "${0##*/}: Invalid setting for TESTING. $TESTING" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Set installation dir.

  local REPO_DIR=$SCRIPT_DIR/..

  if [ -z "${INSTALL_DIR:-}" ] ; then
    # NOTE: this will typiclly not be executed.
    local INSTALL_DIR
    INSTALL_DIR=$REPO_DIR/../install_$(echo $BUILD_TYPE | tr A-Z a-z)
  fi

  #============================================================================
  #---Get mpi stub library if needed.

  local BUILD_DIR=$PWD

  if [ $USE_MPI = OFF ] ; then
    echo "Building mpi-stub ..."
    ln -s ../genomics_gpu/tpls/mpi-stub.tar.gz
    rm -rf mpi-stub
    gunzip <mpi-stub.tar.gz | tar xf -
    pushd mpi-stub
    CC=$CC_serial # NOTE: redefinition!
    make CC=$CC
    popd
    COMET_MPI_COMPILE_OPTS="-I$BUILD_DIR/mpi-stub/include"
    COMET_MPI_LINK_OPTS="-L$BUILD_DIR/mpi-stub/lib -lmpi"
  fi

  #----------------------------------------------------------------------------
  #---Create magma variants.

  if [ $USE_MAGMA = ON ] ; then
    local MAGMA_DIR=$BUILD_DIR/magma_patch
    if [ ! -e $MAGMA_DIR/copy_is_complete ] ; then
      rm -rf $MAGMA_DIR
      echo "Copying magma ..."
      cp -r $REPO_DIR/magma_patch $MAGMA_DIR
      # copy MAGMA source since link will be broken.
      rm $MAGMA_DIR//magma-*.tar.gz
      cp $REPO_DIR/tpls/magma-*.tar.gz $MAGMA_DIR/
      pushd $MAGMA_DIR
      ./create_modified_magmas.sh
      popd
      touch $MAGMA_DIR/copy_is_complete
    fi
  fi

  #----------------------------------------------------------------------------
  #---Compile magma variants.

  if [ $USE_MAGMA = ON ] ; then
    local COMET_MAGMA_COMPILE_OPTS=""
    local COMET_MAGMA_LINK_OPTS=""
    local TAG
    for TAG in minproduct mgemm2 mgemm3 mgemm4 mgemm5 ; do
      local MAGMA_VERSION=magma_$TAG
      local MAGMA_SUBDIR=$MAGMA_DIR/$MAGMA_VERSION
      if [ ! -e $MAGMA_SUBDIR/build_is_complete ] ; then
        pushd $MAGMA_SUBDIR
        ../make_magma.sh
        popd
        if [ -e $MAGMA_SUBDIR/lib/lib${MAGMA_VERSION}.a ] ; then
          touch $MAGMA_SUBDIR/build_is_complete
          get_modules_used_magma > $MAGMA_SUBDIR/modules_used
        else
          exit 1
        fi
      else # build_is_complete
        if [ "$(cksum <$MAGMA_SUBDIR/modules_used)" != \
             "$(cksum < <(get_modules_used_magma))" ] ; then
          echo "${0##*/}: inconsistent modules; please rebuild MAGMA." 1>&2
          exit 1
        fi
      fi # build_is_complete
      COMET_MAGMA_COMPILE_OPTS+=" -I$MAGMA_SUBDIR/include"
      COMET_MAGMA_LINK_OPTS+=" -L$MAGMA_SUBDIR/lib -lmagma_$TAG"
    done # TAG
  fi

  #----------------------------------------------------------------------------
  #---Get unit test harness if needed.

  if [ $TESTING = ON ] ; then
    ln -s ../genomics_gpu/tpls/googletest-release-1.7.0.tar.gz
    # wget -O googletest-release-1.7.0.tar.gz \
    #   https://github.com/google/googletest/archive/release-1.7.0.tar.gz
    gunzip <googletest-release-1.7.0.tar.gz | tar xf -
    local GTEST_DIR=$BUILD_DIR/googletest-release-1.7.0
    mkdir $GTEST_DIR/lib
    $CC_serial -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
      -pthread -c ${GTEST_DIR}/src/gtest-all.cc
    ar -rv $GTEST_DIR/lib/libgtest.a gtest-all.o
    local COMET_TEST_COMPILE_OPTS="-isystem $GTEST_DIR/include -pthread"
    local COMET_TEST_LINK_OPTS="-L$GTEST_DIR/lib -lgtest"
  fi

  #============================================================================
  #---Set variables for cmake.

  local CMAKE_CXX_FLAGS="-DFP_PRECISION_$FP_PRECISION -DADD_"
  CMAKE_CXX_FLAGS+=" -g" # for stack trace
  CMAKE_CXX_FLAGS+=" -Wall -Wno-unused-function -Werror"
  CMAKE_CXX_FLAGS+=" -fno-associative-math"
  CMAKE_CXX_FLAGS+=" -Wno-error=unknown-pragmas"
  CMAKE_CXX_FLAGS+=" -DTEST_PROCS_MAX=$COMET_TEST_PROCS_MAX"
  if [ $USE_GCC = ON ] ; then
      CMAKE_CXX_FLAGS+=" -rdynamic" # for stack trace
  fi

  [[ ${USE_MAGMA:-} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_MAGMA"
  CMAKE_CXX_FLAGS+=" ${COMET_MAGMA_COMPILE_OPTS:-}"

  [[ ${USE_CUDA:-} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_CUDA"
  CMAKE_CXX_FLAGS+=" ${COMET_CUDA_COMPILE_OPTS:-}"

  [[ ${USE_HIP:-} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_HIP"
  CMAKE_CXX_FLAGS+=" ${COMET_HIP_COMPILE_OPTS:-}"

  [[ ${USE_OPENMP:-OFF} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_OPENMP"
  CMAKE_CXX_FLAGS+=" ${COMET_OPENMP_COMPILE_OPTS:-}"

  [[ ${USE_CPUBLAS:-OFF} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_CPUBLAS"
  CMAKE_CXX_FLAGS+=" ${COMET_CPUBLAS_COMPILE_OPTS:-}"

  [[ ${USE_MPI} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_MPI"
  CMAKE_CXX_FLAGS+=" ${COMET_MPI_COMPILE_OPTS:-}"

  [[ ${TESTING:-OFF} = ON ]] && CMAKE_CXX_FLAGS+=" -DTESTING"
  CMAKE_CXX_FLAGS+=" ${COMET_TEST_COMPILE_OPTS:-}"

  [[ ${USE_INT128:-OFF} = ON ]] && CMAKE_CXX_FLAGS+=" -DUSE_INT128"

  CMAKE_CXX_FLAGS+=" ${COMET_EXTRA_COMPILE_OPTS:-}"

  #----------------------------------------------------------------------------
  #---Compiler optimization flags.

  [[ $BUILD_TYPE != "Debug" ]] && CMAKE_CXX_FLAGS+=" -DNDEBUG"
  CMAKE_CXX_FLAGS+=" -O3 -fomit-frame-pointer"

  #---The following change slows performance by 1% on a test case but helps
  #---make results exactly reproducible on varyng number of procs.
  #CMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS -ffast-math"
  CMAKE_CXX_FLAGS+=" -fno-math-errno -ffinite-math-only"
  CMAKE_CXX_FLAGS+=" -fno-signed-zeros -fno-trapping-math -freciprocal-math"
  if [ $USE_GCC = ON ] ; then
    CMAKE_CXX_FLAGS+=" -fno-rounding-math"
    CMAKE_CXX_FLAGS+=" -fno-signaling-nans"
    CMAKE_CXX_FLAGS+=" -fcx-limited-range"
  fi

  CMAKE_CXX_FLAGS+=" -finline-functions"
  if [ $USE_GCC = ON ] ; then
    CMAKE_CXX_FLAGS+=" -finline-limit=1000"
  fi

  #CMAKE_CXX_FLAGS+=" -fstrict-aliasing -fargument-noalias-anything"

  #----------------------------------------------------------------------------
  #---Load flags.

  #local LFLAGS="-v"
  local LFLAGS=""
  LFLAGS+=" ${COMET_MAGMA_LINK_OPTS:-}"
  LFLAGS+=" ${COMET_CUDA_LINK_OPTS:-}"
  LFLAGS+=" ${COMET_HIP_LINK_OPTS:-}"
  LFLAGS+=" ${COMET_CPUBLAS_LINK_OPTS:-}"
  LFLAGS+=" ${COMET_MPI_LINK_OPTS:-}"
  LFLAGS+=" ${COMET_TEST_LINK_OPTS:-}"
  LFLAGS+=" ${COMET_EXTRA_LINK_OPTS:-}"

  #----------------------------------------------------------------------------
  #---Other.

  local CMAKE_EXTRA_OPTIONS=""

#FIX
  if [ "$USE_MPI" = ON ] ; then
    if [ $COMET_PLATFORM = CRAY_XK7 -o $COMET_PLATFORM = EDISON ] ; then
      CMAKE_EXTRA_OPTIONS+=" \
        -DMPI_C_COMPILER="$cc" \
        -DMPI_C_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include \
        -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib \
        -DMPI_CXX_COMPILER="$CC" \
        -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include \
        -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib \
      "
    fi
  fi

#FIX
  if [ ${USE_CUDA:-OFF} = ON ] ; then
    CMAKE_EXTRA_OPTIONS+=" -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
    CMAKE_EXTRA_OPTIONS+=" -DCUDA_HOST_COMPILER:STRING=$(dirname $(which $CC_serial))"
  fi

  #============================================================================
  #---Perform cmake.

  time cmake \
   \
    -DCMAKE_BUILD_TYPE:STRING="$BUILD_TYPE" \
    -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL_DIR" \
   \
    -DCMAKE_C_COMPILER:STRING="$cc" \
    -DCMAKE_CXX_COMPILER:STRING="$CC" \
   \
    -DCMAKE_CXX_FLAGS:STRING="$CMAKE_CXX_FLAGS" \
    -DCMAKE_CXX_FLAGS_DEBUG:STRING="-g -ftrapv" \
    -DCMAKE_CXX_FLAGS_RELEASE:STRING="" \
   \
    -DCMAKE_EXE_LINKER_FLAGS:STRING="$LFLAGS" \
   \
    -DTESTING:BOOL="$TESTING" \
    -DTEST_COMMAND="${COMET_TEST_COMMAND:-}" \
   \
    $CMAKE_EXTRA_OPTIONS \
    -DUSE_MPI:BOOL=${USE_MPI:-OFF} \
    -DUSE_CUDA:BOOL=${USE_CUDA:-OFF} \
    -DUSE_HIP:BOOL=${USE_HIP:-OFF} \
   \
    $REPO_DIR

  ln -s $INSTALL_DIR install_dir
}

#    -DCUDA_HOST_COMPILER:STRING="$(dirname $(which $CC_serial))" \

#==============================================================================

main "$@" 2>&1 | tee out_cmake.txt

#==============================================================================
# cruft.

#  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
#  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
#  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \

  #  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v$DEBUG_FLAG" \
  #  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \

#==============================================================================
