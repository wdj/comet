#!/bin/bash -l
#==============================================================================
#
# Configure script for cmake.
#
# Usage: to configure and build code, perform the following steps:
#
# ...
# cd genomics_gpu
# cd ..
# mkdir build
# cd build
# ../genomics_gpu/scripts/cmake.sh # configure
# ../genomics_gpu/scripts/make.sh # make
#
# Relevant input variables:
#
# REPO_DIR - pathname of cloned CoMet repository
# INSTALL_DIR - installation directory (default - see below)
# NOMPI - build with MPI stub library for single proces use (default OFF)
# BUILD_TYPE - Debug or Release
# FP_PRECISION - SINGLE or DOUBLE precision for floating point operations
#    (default DOUBLE) (does not affect use of tensor cores for CCC)
# TESTING - ON or OFF, to build testing code (default OFF)
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  local IS_CRAY_XK7 # OLCF Titan or Chester
  [[ -n "${CRAYOS_VERSION:-}" ]] && IS_CRAY_XK7="YES" || IS_CRAY_XK7="NO"
  local IS_IBM_AC922 # OLCF Summit or Peak
  [[ -n "${LSF_BINDIR:-}" ]] && IS_IBM_AC922="YES" || IS_IBM_AC922="NO"

  local BUILD_DIR=$PWD

  if [ -z "${REPO_DIR:-}" ] ; then
    local REPO_DIR=$BUILD_DIR/../genomics_gpu
  fi

  #----------------------------------------------------------------------------
  #---Load modules etc.

  if [ $IS_CRAY_XK7 = YES ] ; then
    if [ "$PE_ENV" = "PGI" ] ; then
      module unload PrgEnv-pgi
    fi
    module load PrgEnv-gnu
    module load cudatoolkit
    module load acml
    module load cmake
    local CUDA_INCLUDE_OPTS=$CRAY_CUDATOOLKIT_INCLUDE_OPTS
    local CUDA_POST_LINK_OPTS=$CRAY_CUDATOOLKIT_POST_LINK_OPTS
    local cc=$(which cc)
    local CC=$(which CC)
  elif [ $IS_IBM_AC922 = YES ] ; then
    module -q load gcc/6.4.0
    local CUDA_MODULE=cuda
    module -q load $CUDA_MODULE
    module -q load cmake
    local CUDA_INCLUDE_OPTS="-I$CUDA_DIR/include -I$CUDA_DIR/extras/CUPTI/include -I$CUDA_DIR/extras/Debugger/include"
    local CUDA_POST_LINK_OPTS="-L$CUDA_DIR/targets/ppc64le-linux/lib"
    local cc=$(which mpicc)
    local CC=$(which mpiCC)
  else
    echo "Unknown platform." 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Cleanup of old files.

  rm -rf CMakeCache.txt
  rm -rf CMakeFiles

  #----------------------------------------------------------------------------
  #---Set build type.

  if [ -z "${BUILD_TYPE:-}" ] ; then
    local BUILD_TYPE=Debug #=Release
  fi
  if [ "$BUILD_TYPE" != "Debug" -a "$BUILD_TYPE" != "Release" ] ; then
    echo "Invalid setting for BUILD_TYPE. $BUILD_TYPE" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  # Set whether to build with MPI stub library.

  if [ -z "${NOMPI:-}" ] ; then
    local NOMPI=OFF #=ON
  fi
  if [ "$NOMPI" != "ON" -a "$NOMPI" != "OFF" ] ; then
    echo "Invalid setting for NOMPI. $FP_PRECISION" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Set floating point precision for certain calculations.

  if [ -z "${FP_PRECISION:-}" ] ; then
    local FP_PRECISION=DOUBLE #=SINGLE
  fi
  if [ "$FP_PRECISION" != "SINGLE" -a \
       "$FP_PRECISION" != "DOUBLE" ] ; then
    echo "Invalid setting for FP_PRECISION. $FP_PRECISION" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Set whether to build unit test code.

  if [ -z "${TESTING:-}" ] ; then
    local TESTING=OFF #=ON
  fi
  if [ "$TESTING" != "ON" -a "$TESTING" != "OFF" ] ; then
    echo "Invalid setting for TESTING. $TESTING" 1>&2
    exit 1
  fi

  #----------------------------------------------------------------------------
  #---Set installation dir.

  if [ -z "${INSTALL_DIR:-}" ] ; then
    if [ $BUILD_TYPE = "Debug" ] ; then
      local INSTALL_DIR=$REPO_DIR/../install_debug
    else
      local INSTALL_DIR=$REPO_DIR/../install_release
    fi
  fi

  #============================================================================
  #---Get unit test harness if needed.

  local GTEST_DIR=""

  if [ $TESTING = ON ] ; then
    ln -s ../genomics_gpu/tpls/googletest-release-1.7.0.tar.gz
    # wget -O googletest-release-1.7.0.tar.gz \
    #   https://github.com/google/googletest/archive/release-1.7.0.tar.gz
    gunzip <googletest-release-1.7.0.tar.gz | tar xf -
    GTEST_DIR=$BUILD_DIR/googletest-release-1.7.0
    mkdir $GTEST_DIR/lib
    g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
      -pthread -c ${GTEST_DIR}/src/gtest-all.cc
    ar -rv $GTEST_DIR/lib/libgtest.a gtest-all.o
  fi

  #----------------------------------------------------------------------------
  #---Get mpi stub library if needed.

  if [ $NOMPI = ON ] ; then
    echo "Building mpi-stub ..."
    ln -s ../genomics_gpu/tpls/mpi-stub.tar.gz
    rm -rf mpi-stub
    gunzip <mpi-stub.tar.gz | tar xf -
    pushd mpi-stub
    CC=g++ # NOTE: redefinition!
    make CC=$CC
    popd
  fi

  #----------------------------------------------------------------------------
  #---Create magma variants.

  local MAGMA_DIR=$BUILD_DIR/magma
  if [ ! -e $MAGMA_DIR/copy_is_complete ] ; then
    rm -rf $MAGMA_DIR
    echo "Copying magma ..."
    cp -r $REPO_DIR/magma $MAGMA_DIR
    # copy MAGMA source since link will be broken.
    rm $MAGMA_DIR//magma-*.tar.gz
    cp $REPO_DIR/tpls/magma-*.tar.gz $MAGMA_DIR/
    pushd $MAGMA_DIR
    ./create_modified_magmas.sh
    popd
    touch $MAGMA_DIR/copy_is_complete
  fi

  #----------------------------------------------------------------------------
  #---Compile magma variants.

  local tag
  for tag in minproduct tally4 tally3 tally2 ; do
    local magma_version=magma_$tag
    local magma_subdir=$MAGMA_DIR/$magma_version
    if [ ! -e $magma_subdir/build_is_complete ] ; then
      pushd $magma_subdir
      ../make_magma.sh
      popd
      if [ -e $magma_subdir/lib/lib${magma_version}.a ] ; then
        touch $magma_subdir/build_is_complete
      else
        exit 1
      fi
    fi
  done

  #============================================================================
  #---Set variables for cmake.

  local C_CXX_FLAGS
  C_CXX_FLAGS="-DFP_PRECISION_$FP_PRECISION -DADD_"
  C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_minproduct/include"
  C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_tally4/include"
  C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_tally3/include"
  C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_tally2/include"
  C_CXX_FLAGS="$C_CXX_FLAGS $CUDA_INCLUDE_OPTS"
  C_CXX_FLAGS="$C_CXX_FLAGS -g -rdynamic" # for stack trace
  C_CXX_FLAGS="$C_CXX_FLAGS -Wall -Wno-unused-function -Werror"
  C_CXX_FLAGS="$C_CXX_FLAGS -fno-associative-math"
  C_CXX_FLAGS="$C_CXX_FLAGS -fopenmp"
  C_CXX_FLAGS="$C_CXX_FLAGS -Wno-error=unknown-pragmas"
  C_CXX_FLAGS="$C_CXX_FLAGS -DTEST_PROCS_MAX=64"
  #C_CXX_FLAGS="$C_CXX_FLAGS -std=c++11"
  #C_CXX_FLAGS="$C_CXX_FLAGS -Wconversion"

  #----------------------------------------------------------------------------

  local C_FLAGS_OPT="-DNDEBUG -O3 -fomit-frame-pointer"

  #---The following change slows performance by 1% on a test case but helps
  #---make results exactly reproducible on varyng number of procs.
  #C_FLAGS_OPT="$C_FLAGS_OPT -ffast-math"
  C_FLAGS_OPT="$C_FLAGS_OPT -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fcx-limited-range"
  C_FLAGS_OPT="$C_FLAGS_OPT -fno-signed-zeros -fno-trapping-math -freciprocal-math"

  C_FLAGS_OPT="$C_FLAGS_OPT -finline-functions -finline-limit=1000"
  if [ $IS_CRAY_XK7 = YES ] ; then
    C_FLAGS_OPT="$C_FLAGS_OPT -march=bdver1"
  fi
  if [ $IS_IBM_AC922 = YES ] ; then
    C_FLAGS_OPT="$C_FLAGS_OPT -mcpu=power9 -mtune=power9 -mcmodel=large -m64"
  fi
  #C_FLAGS_OPT="$C_FLAGS_OPT -fstrict-aliasing -fargument-noalias-anything"

  local C_FLAGS_RELEASE=""
  #C_FLAGS_RELEASE="$C_FLAGS_OPT"
  local C_CXX_FLAGS="$C_CXX_FLAGS $C_FLAGS_OPT"

  #----------------------------------------------------------------------------

  local LFLAGS=
  LFLAGS="-L$MAGMA_DIR/magma_minproduct/lib -lmagma_minproduct"
  LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally4/lib -lmagma_tally4"
  LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally3/lib -lmagma_tally3"
  LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally2/lib -lmagma_tally2"
  LFLAGS="$LFLAGS $CUDA_POST_LINK_OPTS -lcublas -lcudart"
  if [ $IS_CRAY_XK7 = YES ] ; then
    LFLAGS="$LFLAGS -Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib"
    LFLAGS="$LFLAGS -Wl,-rpath=/opt/acml/5.3.1/gfortran64_mp/lib"
  fi
  if [ $IS_IBM_AC922 = YES ] ; then
    LFLAGS="$LFLAGS -Wl,-rpath=$CUDA_DIR/lib64"
  fi

  #----------------------------------------------------------------------------

  if [ $IS_CRAY_XK7 = YES ] ; then
    local TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64"
    C_CXX_FLAGS="$C_CXX_FLAGS -DHAVE_INT128"
  fi
  if [ $IS_IBM_AC922 = YES ] ; then
    local TEST_COMMAND="module load $CUDA_MODULE ; env OMP_NUM_THREADS=16 jsrun -n 2 -r 1 -c 32 -g 6 -a 32 -X 1"
    C_CXX_FLAGS="$C_CXX_FLAGS -DUSE_TC"
    C_CXX_FLAGS="$C_CXX_FLAGS -DHAVE_INT128"
  fi

  if [ "$NOMPI" = ON ] ; then
    C_CXX_FLAGS="$C_CXX_FLAGS -DNOMPI -I$BUILD_DIR/mpi-stub/include"
    LFLAGS="$LFLAGS -L$PWD/mpi-stub/lib -lmpi"
  fi

  #  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v$DEBUG_FLAG" \
  #  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \

  local CMAKE_MPI_OPTIONS=""
  if [ "$NOMPI" = OFF ] ; then
    if [ $IS_CRAY_XK7 = YES ] ; then
      CMAKE_MPI_OPTIONS="\
        -DMPI_C_COMPILER="$cc" \
        -DMPI_CXX_COMPILER="$CC" \
        -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include \
        -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib \
      "
    fi
  else
    local CMAKE_MPI_OPTIONS="-DNOMPI=ON"
  fi

  #============================================================================
  #---Perform cmake.

  time cmake \
   \
    -DCMAKE_BUILD_TYPE:STRING="$BUILD_TYPE" \
    -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL_DIR" \
   \
    $CMAKE_MPI_OPTIONS \
   \
    -DCMAKE_C_COMPILER:STRING="$cc" \
   \
    -DCMAKE_CXX_COMPILER:STRING="$CC" \
   \
    -DC_AND_CXX_FLAGS:STRING="$C_CXX_FLAGS" \
   \
    -DCMAKE_C_FLAGS:STRING="-std=c99 -pedantic" \
    -DCMAKE_C_FLAGS_DEBUG:STRING="-g -ftrapv" \
    -DCMAKE_C_FLAGS_RELEASE:STRING="$C_FLAGS_RELEASE" \
   \
    -DCMAKE_CXX_FLAGS:STRING="" \
    -DCMAKE_CXX_FLAGS_DEBUG:STRING="-g -ftrapv" \
    -DCMAKE_CXX_FLAGS_RELEASE:STRING="$C_FLAGS_RELEASE" \
   \
    -DCMAKE_EXE_LINKER_FLAGS:STRING="$LFLAGS" \
   \
    -DTESTING:BOOL=$TESTING \
    -DGTEST_DIR:STRING=$GTEST_DIR \
    -DTEST_COMMAND="$TEST_COMMAND" \
   \
    -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \
   \
    $REPO_DIR

	ln -s $INSTALL_DIR install_dir
}

#  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
#  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
#  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \

#==============================================================================

main "$@"

#==============================================================================
