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

#---Modules.

if [ "$PE_ENV" = "PGI" ] ; then
  module swap PrgEnv-pgi PrgEnv-gnu
fi
module load cudatoolkit
module load acml
module load cmake
#module load cmake3/3.2.3

#---Cleanup.

rm -rf CMakeCache.txt
rm -rf CMakeFiles

#---Main project dir for cloned repo.

if [ "$PROJECT_DIR" = "" ] ; then
  PROJECT_DIR=${PWD}/../genomics_gpu
fi

#---Build type.

if [ "$BUILD_TYPE" = "" ] ; then
  BUILD_TYPE=Debug
  #BUILD_TYPE=Release
fi

#---Installation dir.

if [ "$INSTALL_DIR" = "" ] ; then
  if [ "$BUILD_TYPE" = "Debug" ] ; then
    INSTALL_DIR=${PROJECT_DIR}/../install_debug
  else
    INSTALL_DIR=${PROJECT_DIR}/../install_release
  fi
fi

#---Floating point precision for calculations.

if [ "$FP_PRECISION" = "" ] ; then
  #FP_PRECISION=FP_PRECISION_SINGLE
  FP_PRECISION=FP_PRECISION_DOUBLE
fi

#---Whether to build testing code.

if [ "$TESTING" = "" ] ; then
  TESTING=OFF
  #TESTING=ON
fi

#==============================================================================

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
  -DC_AND_CXX_FLAGS:STRING="-D$FP_PRECISION -I$PROJECT_DIR/magma/magma_minproduct-1.6.2/include -DADD_ $CRAY_CUDATOOLKIT_INCLUDE_OPTS" \
 \
  -DCMAKE_C_FLAGS:STRING="-std=c99 -Wall -Wno-unused-function -pedantic -Werror" \
  -DCMAKE_C_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="-DNDEBUG -O3 -ffast-math -fargument-noalias-anything -fstrict-aliasing -finline-functions -finline-limit=1000 -fomit-frame-pointer" \
 \
  -DCMAKE_CXX_FLAGS:STRING="-Wall -Wno-unused-function -Werror" \
  -DCMAKE_CXX_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_CXX_FLAGS_RELEASE:STRING="-DNDEBUG -O3 -ffast-math -fargument-noalias-anything -fstrict-aliasing -finline-functions -finline-limit=1000 -fomit-frame-pointer" \
 \
  -DCMAKE_EXE_LINKER_FLAGS:STRING="-L$PROJECT_DIR/magma/magma_minproduct-1.6.2/lib -lmagma_minproduct $CRAY_CUDATOOLKIT_POST_LINK_OPTS -lcublas" \
 \
  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \
  -DTESTING:BOOL=$TESTING \
  -DGTEST_DIR:STRING=$GTEST_DIR \
 \
  $PROJECT_DIR

#==============================================================================
