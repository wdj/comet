#!/bin/bash -l
#==============================================================================
#
# Configure script for cmake.
#
# Usage: to configure and build code:
#
# ...
# cd genomics_gpu
# cd ..
# mkdir build
# cd build
# ../genomics_gpu/scripts/cmake.sh
# ../genomics_gpu/scripts/make.sh
#
#==============================================================================

#------------------------------------------------------------------------------
#---Load modules.

if [ "$PE_ENV" = "PGI" ] ; then
  module swap PrgEnv-pgi PrgEnv-gnu
fi
module load cudatoolkit
module load acml
module load cmake
#module load cmake3/3.2.3

#------------------------------------------------------------------------------
#---Cleanup.

rm -rf CMakeCache.txt
rm -rf CMakeFiles

#------------------------------------------------------------------------------
#---Main project dir for cloned repo.

if [ "$PROJECT_DIR" = "" ] ; then
  PROJECT_DIR=${PWD}/../genomics_gpu
fi

#------------------------------------------------------------------------------
#---Set build type.

if [ "$BUILD_TYPE" = "" ] ; then
  BUILD_TYPE=Debug
  #BUILD_TYPE=Release
fi
if [ "$BUILD_TYPE" != "Debug" -a "$BUILD_TYPE" != "Release" ] ; then
  echo "Invalid setting for BUILD_TYPE. $BUILD_TYPE" 1>&2
  exit 1
fi

#------------------------------------------------------------------------------
#---Set installation dir.

if [ "$INSTALL_DIR" = "" ] ; then
  if [ "$BUILD_TYPE" = "Debug" ] ; then
    INSTALL_DIR=${PROJECT_DIR}/../install_debug
  else
    INSTALL_DIR=${PROJECT_DIR}/../install_release
  fi
fi

#------------------------------------------------------------------------------
#---Set floating point precision for certain calculations.

if [ "$FP_PRECISION" = "" ] ; then
  #FP_PRECISION=SINGLE
  FP_PRECISION=DOUBLE
fi
if [ "$FP_PRECISION" != "SINGLE" -a \
     "$FP_PRECISION" != "DOUBLE" ] ; then
  echo "Invalid setting for FP_PRECISION. $FP_PRECISION" 1>&2
  exit 1
fi

#------------------------------------------------------------------------------
#---Set whether to build unit test code.

if [ "$TESTING" = "" ] ; then
  TESTING=OFF
  #TESTING=ON
fi
if [ "$TESTING" != "ON" -a "$TESTING" != "OFF" ] ; then
  echo "Invalid setting for TESTING. $TESTING" 1>&2
  exit 1
fi

#==============================================================================
#---Get unit test harness if needed.

GTEST_DIR=""

if [ "$TESTING" = ON ] ; then
  wget -O googletest-release-1.7.0.tar.gz \
    https://github.com/google/googletest/archive/release-1.7.0.tar.gz
  gunzip <googletest-release-1.7.0.tar.gz | tar xf -
  GTEST_DIR=$(pwd)/googletest-release-1.7.0
  mkdir $GTEST_DIR/lib
  g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
    -pthread -c ${GTEST_DIR}/src/gtest-all.cc
  ar -rv $GTEST_DIR/lib/libgtest.a gtest-all.o
fi

#==============================================================================
#---Build magma variants.

MAGMA_DIR=$PWD/magma
if [ ! -e $MAGMA_DIR/copy_is_complete ] ; then
  cp -rp $PROJECT_DIR/magma $MAGMA_DIR
  touch $MAGMA_DIR/copy_is_complete
fi

for magma_version in magma_minproduct magma_tally3 magma_tally4 ; do
  MAGMA_SUBDIR=$MAGMA_DIR/$magma_version
  if [ ! -e $MAGMA_SUBDIR/build_is_complete ] ; then
    pushd $MAGMA_SUBDIR
    ../make.sh
    popd
    touch $MAGMA_SUBDIR/build_is_complete
  fi
done

#==============================================================================
#---Perform cmake.

C_CXX_FLAGS="-DFP_PRECISION_$FP_PRECISION -DADD_"
C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_minproduct/include"
C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_tally4/include"
C_CXX_FLAGS="$C_CXX_FLAGS -I$MAGMA_DIR/magma_tally3/include"
C_CXX_FLAGS="$C_CXX_FLAGS $CRAY_CUDATOOLKIT_INCLUDE_OPTS"

#----------

C_FLAGS_RELEASE="-DNDEBUG -O3 -fomit-frame-pointer"

#---The following change slows performance by 1% on a test case but helps
#---make results exactly reproducible on vryng number of procs.
#C_FLAGS_RELEASE="$C_FLAGS_RELEASE -ffast-math"
C_FLAGS_RELEASE="$C_FLAGS_RELEASE -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fcx-limited-range"
C_FLAGS_RELEASE="$C_FLAGS_RELEASE -fno-signed-zeros -fno-trapping-math -freciprocal-math"

C_FLAGS_RELEASE="$C_FLAGS_RELEASE -finline-functions -finline-limit=1000"
#C_FLAGS_RELEASE="$C_FLAGS_RELEASE -fstrict-aliasing -fargument-noalias-anything"

#----------

LFLAGS="-L$MAGMA_DIR/magma_minproduct/lib -lmagma_minproduct"
LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally4/lib -lmagma_tally4"
LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally3/lib -lmagma_tally3"
LFLAGS="$LFLAGS $CRAY_CUDATOOLKIT_POST_LINK_OPTS -lcublas"
LFLAGS="$LFLAGS -Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib"

#------------------------------------------------------------------------------

time cmake \
 \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD_TYPE" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL_DIR" \
 \
  -DCMAKE_C_COMPILER:STRING="$(which cc)" \
  -DMPI_C_COMPILER="$(which cc)" \
 \
  -DCMAKE_CXX_COMPILER:STRING="$(which CC)" \
  -DMPI_CXX_COMPILER="$(which CC)" \
  -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include \
  -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib \
 \
  -DC_AND_CXX_FLAGS:STRING="$C_CXX_FLAGS" \
 \
  -DCMAKE_C_FLAGS:STRING="-std=c99 -Wall -Wno-unused-function -pedantic -Werror -fno-associative-math" \
  -DCMAKE_C_FLAGS_DEBUG:STRING="-g -ftrapv" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="$C_FLAGS_RELEASE" \
 \
  -DCMAKE_CXX_FLAGS:STRING="-Wall -Wno-unused-function -Werror -fno-associative-math" \
  -DCMAKE_CXX_FLAGS_DEBUG:STRING="-g -ftrapv" \
  -DCMAKE_CXX_FLAGS_RELEASE:STRING="$C_FLAGS_RELEASE" \
 \
  -DCMAKE_EXE_LINKER_FLAGS:STRING="$LFLAGS" \
 \
  -DTESTING:BOOL=$TESTING \
  -DGTEST_DIR:STRING=$GTEST_DIR \
 \
  $PROJECT_DIR

ln -s $INSTALL_DIR install_dir

#  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
#  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
#  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \

#==============================================================================
