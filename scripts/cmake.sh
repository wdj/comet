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

BUILD_DIR=$PWD

if [ -n "${CRAYOS_VERSION:-}" ] ; then
  IS_CRAY_XK7="YES"
else
  IS_CRAY_XK7="NO"
fi

#------------------------------------------------------------------------------
#---Load modules.

if [ $IS_CRAY_XK7 = YES ] ; then
  # For Titan or Chester
  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  #module swap gcc gcc/4.9.3
  module load cudatoolkit
  #module load cudatoolkit/7.5.18-1.0502.10743.2.1
  module load acml
  module load cmake
  CUDA_INCLUDE_OPTS=$CRAY_CUDATOOLKIT_INCLUDE_OPTS
  CUDA_POST_LINK_OPTS=$CRAY_CUDATOOLKIT_POST_LINK_OPTS
  cc=$(which cc)
  CC=$(which CC)
else #---IBM
  # For Summit or Peak
  #module load spectrum-mpi/10.2.0.0-20180508 #FIX
  module load gcc/6.4.0
  #CUDA_MODULE=cuda/9.1.85
  CUDA_MODULE=cuda
  module load $CUDA_MODULE
  module load cmake
  CUDA_INCLUDE_OPTS="-I$CUDA_DIR/include -I$CUDA_DIR/extras/CUPTI/include -I$CUDA_DIR/extras/Debugger/include"
  CUDA_POST_LINK_OPTS="-L$CUDA_DIR/targets/ppc64le-linux/lib"
  cc=$(which mpicc)
  CC=$(which mpiCC)
fi

#------------------------------------------------------------------------------
#---Cleanup.

rm -rf CMakeCache.txt
rm -rf CMakeFiles

#------------------------------------------------------------------------------
#---Main project dir for cloned repo.

if [ -z "$REPO_DIR" ] ; then
  REPO_DIR=$BUILD_DIR/../genomics_gpu
fi

#------------------------------------------------------------------------------
#---Set build type.

if [ -z "$BUILD_TYPE" ] ; then
  BUILD_TYPE=Debug
  #BUILD_TYPE=Release
fi
if [ "$BUILD_TYPE" != "Debug" -a "$BUILD_TYPE" != "Release" ] ; then
  echo "Invalid setting for BUILD_TYPE. $BUILD_TYPE" 1>&2
  exit 1
fi

#------------------------------------------------------------------------------
# Set whether to build with MPI stub library.

if [ -z "$NOMPI" ] ; then
  NOMPI=OFF
fi
if [ "$NOMPI" != "ON" -a "$NOMPI" != "OFF" ] ; then
  echo "Invalid setting for NOMPI. $FP_PRECISION" 1>&2
  exit 1
fi

#------------------------------------------------------------------------------
#---Set floating point precision for certain calculations.

if [ -z "$FP_PRECISION" ] ; then
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

if [ -z "$TESTING" ] ; then
  TESTING=OFF
  #TESTING=ON
fi
if [ "$TESTING" != "ON" -a "$TESTING" != "OFF" ] ; then
  echo "Invalid setting for TESTING. $TESTING" 1>&2
  exit 1
fi

#------------------------------------------------------------------------------
#---Set installation dir.

if [ -z "$INSTALL_DIR" ] ; then
  if [ "$BUILD_TYPE" = "Debug" ] ; then
    INSTALL_DIR=$REPO_DIR/../install_debug
  else
    INSTALL_DIR=$REPO_DIR/../install_release
  fi
fi

#==============================================================================
#---Get unit test harness if needed.

GTEST_DIR=""

if [ "$TESTING" = ON ] ; then
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

#==============================================================================
#---Get mpi stub library if needed.

if [ "$NOMPI" = ON ] ; then
  echo "Building mpi-stub ..."
  ln -s ../genomics_gpu/tpls/mpi-stub.tar.gz
  rm -rf mpi-stub
  gunzip <mpi-stub.tar.gz | tar xf -
  pushd mpi-stub
  CC=g++ # NOTE: redefinitin!
  make CC=$CC
  popd
fi

#==============================================================================
#---Build magma variants.

MAGMA_DIR=$BUILD_DIR/magma
if [ ! -e $MAGMA_DIR/copy_is_complete ] ; then
  rm -rf $MAGMA_DIR
  echo "Copying magma ..."
  cp -r $REPO_DIR/magma $MAGMA_DIR
  rm $MAGMA_DIR//magma-*.tar.gz
  cp $REPO_DIR/tpls/magma-*.tar.gz $MAGMA_DIR/
  pushd $MAGMA_DIR
  ./create_modified_magmas.sh
  popd
  touch $MAGMA_DIR/copy_is_complete
fi

for magma_version in magma_minproduct magma_tally4 magma_tally3 \
                     magma_tally2 ; do
  MAGMA_SUBDIR=$MAGMA_DIR/$magma_version
  if [ ! -e $MAGMA_SUBDIR/build_is_complete ] ; then
    pushd $MAGMA_SUBDIR
    ../make_magma.sh
    popd
    if [ -e $MAGMA_SUBDIR/lib/lib${magma_version}.a ] ; then
      touch $MAGMA_SUBDIR/build_is_complete
    else
      exit 1
    fi
  fi
done

#==============================================================================
#---Perform cmake.

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

#----------

C_FLAGS_OPT="-DNDEBUG -O3 -fomit-frame-pointer"

#---The following change slows performance by 1% on a test case but helps
#---make results exactly reproducible on varyng number of procs.
#C_FLAGS_OPT="$C_FLAGS_OPT -ffast-math"
C_FLAGS_OPT="$C_FLAGS_OPT -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fcx-limited-range"
C_FLAGS_OPT="$C_FLAGS_OPT -fno-signed-zeros -fno-trapping-math -freciprocal-math"

C_FLAGS_OPT="$C_FLAGS_OPT -finline-functions -finline-limit=1000"
if [ $IS_CRAY_XK7 = YES ] ; then
  C_FLAGS_OPT="$C_FLAGS_OPT -march=bdver1"
else
  C_FLAGS_OPT="$C_FLAGS_OPT -mcpu=power9 -mtune=power9 -mcmodel=large -m64"
fi
#C_FLAGS_OPT="$C_FLAGS_OPT -fstrict-aliasing -fargument-noalias-anything"

C_FLAGS_RELEASE=""
#C_FLAGS_RELEASE="$C_FLAGS_OPT"
C_CXX_FLAGS="$C_CXX_FLAGS $C_FLAGS_OPT"

#----------

LFLAGS="-L$MAGMA_DIR/magma_minproduct/lib -lmagma_minproduct"
LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally4/lib -lmagma_tally4"
LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally3/lib -lmagma_tally3"
LFLAGS="$LFLAGS -L$MAGMA_DIR/magma_tally2/lib -lmagma_tally2"
LFLAGS="$LFLAGS $CUDA_POST_LINK_OPTS -lcublas -lcudart"
if [ $IS_CRAY_XK7 = YES ] ; then
  LFLAGS="$LFLAGS -Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib"
  LFLAGS="$LFLAGS -Wl,-rpath=/opt/acml/5.3.1/gfortran64_mp/lib"
else
  LFLAGS="$LFLAGS -Wl,-rpath=$CUDA_DIR/lib64"
fi

if [ $IS_CRAY_XK7 = YES ] ; then
  TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64"
  C_CXX_FLAGS="$C_CXX_FLAGS -DHAVE_INT128"
else
  TEST_COMMAND="module load $CUDA_MODULE ; env OMP_NUM_THREADS=16 jsrun -n 2 -r 1 -c 32 -g 6 -a 32 -X 1"
  C_CXX_FLAGS="$C_CXX_FLAGS -DUSE_TC"
  C_CXX_FLAGS="$C_CXX_FLAGS -DHAVE_INT128"
fi

if [ "$NOMPI" = ON ] ; then
  C_CXX_FLAGS="$C_CXX_FLAGS -DNOMPI -I$BUILD_DIR/mpi-stub/include"
  LFLAGS="$LFLAGS -L$PWD/mpi-stub/lib -lmpi"
fi

#  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v$DEBUG_FLAG" \
#  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \

#------------------------------------------------------------------------------

if [ "$NOMPI" = OFF ] ; then
  CMAKE_MPI_OPTIONS="\
    -DMPI_C_COMPILER="$cc" \
    -DMPI_CXX_COMPILER="$CC" \
    -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include \
    -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib \
  "
else
  CMAKE_MPI_OPTIONS="-DNOMPI=ON"
fi

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

#  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
#  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
#  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \

#==============================================================================
