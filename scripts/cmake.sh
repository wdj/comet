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

  local USE_CUDA=ON
  local USE_MAGMA=ON
  local USE_HIP=OFF
  local USE_CPUBLAS=OFF
  local CC_serial=g++

  if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
    true # skip
  #--------------------
  elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    local CUDA_INCLUDE_OPTS=$CRAY_CUDATOOLKIT_INCLUDE_OPTS
    local CUDA_POST_LINK_OPTS=$CRAY_CUDATOOLKIT_POST_LINK_OPTS
    local cc=$(which cc)
    local CC=$(which CC)
  #--------------------
  elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
    local CUDA_ROOT=$OLCF_CUDA_ROOT
    local CUDA_INCLUDE_OPTS="-I$CUDA_ROOT/include"
    CUDA_INCLUDE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
    CUDA_INCLUDE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
    local CUDA_POST_LINK_OPTS="-L$OLCF_CUDA_ROOT/targets/ppc64le-linux/lib"
    local cc=$(which mpicc)
    local CC=$(which mpiCC)
    USE_CPUBLAS=ON
  #--------------------
  elif [ $COMET_PLATFORM = DGX2 ] ; then
    local CUDA_ROOT="$HOME/cuda"
    local CUDA_INCLUDE_OPTS="-I$CUDA_ROOT/include"
    CUDA_INCLUDE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
    CUDA_INCLUDE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
    local CUDA_POST_LINK_OPTS="-L$CUDA_ROOT/lib64"
    local cc=$HOME/.linuxbrew/bin/gcc-6
    local CC=$HOME/.linuxbrew/bin/g++-6
  #--------------------
  elif [ $COMET_PLATFORM = GPUSYS2 ] ; then
    export CUDA_ROOT=/usr/local/cuda-10.1
    local CUDA_INCLUDE_OPTS="-I$CUDA_ROOT/include"
    CUDA_INCLUDE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
    CUDA_INCLUDE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
    local CUDA_POST_LINK_OPTS="-L$CUDA_ROOT/lib64"
    local cc=$(spack location --install-dir gcc)/bin/gcc
    local CC=$(spack location --install-dir gcc)/bin/g++
    local CC_serial=$CC
  #--------------------
  elif [ $COMET_PLATFORM = EDISON ] ; then
    USE_CUDA=OFF
    USE_MAGMA=OFF
    local CUDA_INCLUDE_OPTS=""
    local CUDA_POST_LINK_OPTS=""
    local cc=$(which cc)
    local CC=$(which CC)
  #--------------------
  elif [ $COMET_PLATFORM = LYRA ] ; then
    USE_CUDA=OFF
    USE_MAGMA=OFF
    USE_HIP=ON
    local CUDA_INCLUDE_OPTS=""
    local CUDA_POST_LINK_OPTS=""
    local cc=$(which gcc)
    local CC=$(which g++)
    CC_serial=hipcc
    if [ "${BLIS_PATH:-}" != "" ] ; then
      USE_CPUBLAS=ON
    fi
  #--------------------
  elif [ $COMET_PLATFORM = AMDINTERNAL ] ; then
    USE_CUDA=OFF
    USE_MAGMA=OFF
    USE_HIP=ON
    local CUDA_INCLUDE_OPTS=""
    local CUDA_POST_LINK_OPTS=""
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

  if [ $USE_MPI = OFF ] ; then
    echo "Building mpi-stub ..."
    ln -s ../genomics_gpu/tpls/mpi-stub.tar.gz
    rm -rf mpi-stub
    gunzip <mpi-stub.tar.gz | tar xf -
    pushd mpi-stub
    CC=$CC_serial # NOTE: redefinition!
    make CC=$CC
    popd
  fi

  #----------------------------------------------------------------------------
  #---Create magma variants.

  local BUILD_DIR=$PWD

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
    done # TAG
  fi

  #----------------------------------------------------------------------------
  #---Get unit test harness if needed.

  local GTEST_DIR=""
  local GTEST_LD_FLAGS=""

  if [ $TESTING = ON ] ; then
    ln -s ../genomics_gpu/tpls/googletest-release-1.7.0.tar.gz
    # wget -O googletest-release-1.7.0.tar.gz \
    #   https://github.com/google/googletest/archive/release-1.7.0.tar.gz
    gunzip <googletest-release-1.7.0.tar.gz | tar xf -
    GTEST_DIR=$BUILD_DIR/googletest-release-1.7.0
    mkdir $GTEST_DIR/lib
    $CC_serial -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
      -pthread -c ${GTEST_DIR}/src/gtest-all.cc
    ar -rv $GTEST_DIR/lib/libgtest.a gtest-all.o
    GTEST_LD_FLAGS="-L$GTEST_DIR/lib -lgtest"
  fi

  #============================================================================
  #---Set variables for cmake.

  local CMAKE_CXX_FLAGS="-DFP_PRECISION_$FP_PRECISION -DADD_"
  CMAKE_CXX_FLAGS+=" -g" # for stack trace
  if [ $USE_HIP = OFF ] ; then
    CMAKE_CXX_FLAGS+=" -rdynamic" # for stack trace
  fi
  CMAKE_CXX_FLAGS+=" -Wall -Wno-unused-function -Werror"
  CMAKE_CXX_FLAGS+=" -fno-associative-math"
  if [ $USE_HIP = OFF ] ; then
    CMAKE_CXX_FLAGS+=" -fopenmp" # TODO: put this back in later
  fi
  CMAKE_CXX_FLAGS+=" -Wno-error=unknown-pragmas"
  CMAKE_CXX_FLAGS+=" -DTEST_PROCS_MAX=64"
  #CMAKE_CXX_FLAGS+=" -std=c++11 -Wconversion"

  if [ $USE_CUDA = ON ] ; then
    CMAKE_CXX_FLAGS+=" -DUSE_CUDA $CUDA_INCLUDE_OPTS"
  fi

  if [ $USE_MAGMA = ON ] ; then
    CMAKE_CXX_FLAGS+=" -DUSE_MAGMA"
    local TAG
    for TAG in minproduct mgemm2 mgemm3 mgemm4 mgemm5 ; do
      CMAKE_CXX_FLAGS+=" -I$MAGMA_DIR/magma_$TAG/include"
    done
  fi

  if [ $USE_HIP = ON ] ; then
    CMAKE_CXX_FLAGS+=" -DUSE_HIP"
    CMAKE_CXX_FLAGS+=" -I$ROCBLAS_PATH/include"
    CMAKE_CXX_FLAGS+=" -I$ROCM_PATH/include"
    CMAKE_CXX_FLAGS+=" -I$HIP_PATH/include/hip"
    #CMAKE_CXX_FLAGS+=" -D__HIP_PLATFORM_HCC__"
  fi

  if [ $USE_CPUBLAS = ON ] ; then
    CMAKE_CXX_FLAGS+=" -DUSE_CPUBLAS"
    if [ $COMET_PLATFORM = IBM_AC922 ] ; then
      CMAKE_CXX_FLAGS+=" -I$OLCF_ESSL_ROOT/include"
    fi
    if [ $COMET_PLATFORM = LYRA ] ; then
      CMAKE_CXX_FLAGS+=" -I$BLIS_PATH/include/zen"
    fi
  fi

  if [ $TESTING = ON ] ; then
    CMAKE_CXX_FLAGS+=" -isystem $GTEST_DIR/include -pthread -DTESTING"
  fi

  #----------------------------------------------------------------------------
  #---Compiler optimization flags.

  local CXX_FLAGS_OPT="-O3 -fomit-frame-pointer"
  [[ $BUILD_TYPE != "Debug" ]] && CXX_FLAGS_OPT+=" -DNDEBUG"

  #---The following change slows performance by 1% on a test case but helps
  #---make results exactly reproducible on varyng number of procs.
  #CXX_FLAGS_OPT="$CXX_FLAGS_OPT -ffast-math"
  CXX_FLAGS_OPT+=" -fno-math-errno -ffinite-math-only"
  if [ $USE_HIP = OFF ] ; then
    CXX_FLAGS_OPT+=" -fno-rounding-math"
    CXX_FLAGS_OPT+=" -fno-signaling-nans"
    CXX_FLAGS_OPT+=" -fcx-limited-range"
  fi
  CXX_FLAGS_OPT+=" -fno-signed-zeros -fno-trapping-math -freciprocal-math"

  CXX_FLAGS_OPT+=" -finline-functions"
  if [ $USE_HIP = OFF ] ; then
    CXX_FLAGS_OPT+=" -finline-limit=1000"
  fi
  if [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    CXX_FLAGS_OPT+=" -march=bdver1"
  fi
  if [ $COMET_PLATFORM = IBM_AC922 ] ; then
    CXX_FLAGS_OPT+=" -mcpu=power9 -mtune=power9 -mcmodel=large -m64"
  fi
  #CXX_FLAGS_OPT+=" -fstrict-aliasing -fargument-noalias-anything"

  CMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS $CXX_FLAGS_OPT"

  #----------------------------------------------------------------------------
  #---Load flags.

  local LFLAGS=""
  if [ $USE_MAGMA = ON ] ; then
    local TAG
    for TAG in minproduct mgemm2 mgemm3 mgemm4 mgemm5 ; do
      LFLAGS+="-L$MAGMA_DIR/magma_$TAG/lib -lmagma_$TAG "
    done
  fi

  if [ $USE_CUDA = ON ] ; then
    LFLAGS+=" $CUDA_POST_LINK_OPTS -lcublas -lcudart"
  fi

  if [ $USE_HIP = ON ] ; then
    LFLAGS+=" -L$ROCBLAS_PATH/lib -lrocblas"
    LFLAGS+=" -L$ROCM_PATH/lib -lhip_hcc"
    if [ $USE_CPUBLAS = ON ] ; then
      LFLAGS+=" -L$BLIS_PATH/lib/zen"
      LFLAGS+=" -Wl,-rpath,$BLIS_PATH/lib/zen -lblis"
    fi
  fi

  if [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    LFLAGS+=" -Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib"
    LFLAGS+=" -Wl,-rpath=/opt/acml/5.3.1/gfortran64_mp/lib"
  fi

  if [ $COMET_PLATFORM = IBM_AC922 ] ; then
    LFLAGS+=" -Wl,-rpath=$OLCF_CUDA_ROOT/lib64"
    if [ $USE_CPUBLAS = ON ] ; then
      local XLF_DIR=$(module load xl 2>/dev/null ; echo $OLCF_XLF_ROOT)/lib
      local XLF_DIR2=$(module load xl 2>/dev/null ; echo $OLCF_XL_ROOT)/lib
      LFLAGS+=" -L$OLCF_ESSL_ROOT/lib64"
      LFLAGS+=" -Wl,-rpath,$OLCF_ESSL_ROOT/lib64 -lessl"
      LFLAGS+=" -L$XLF_DIR -Wl,-rpath,$XLF_DIR2 -lxlf90_r"
      LFLAGS+=" -lxl -lxlfmath"
      LFLAGS+=" -Wl,-rpath,$OLCF_GCC_ROOT/lib64"
    fi
  fi

  if [ $COMET_PLATFORM = DGX2 -o $COMET_PLATFORM = GPUSYS2 ] ; then
    CMAKE_CXX_FLAGS+=" -std=gnu++11"
    LFLAGS+=" -Wl,-rpath=$CUDA_ROOT/lib64"
  fi

  if [ $COMET_PLATFORM = EDISON ] ; then
    CMAKE_CXX_FLAGS+=" -std=gnu++11"
  fi

  if [ "$USE_MPI" = ON ] ; then
    CMAKE_CXX_FLAGS+=" -DUSE_MPI"
  else
    CMAKE_CXX_FLAGS+=" -I$BUILD_DIR/mpi-stub/include"
    LFLAGS+=" -L$PWD/mpi-stub/lib -lmpi"
  fi

  #----------------------------------------------------------------------------
  #---Test flags.

  local TEST_COMMAND=""

  if [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64"
    CMAKE_CXX_FLAGS+=" -DHAVE_INT128"
  fi

  if [ $COMET_PLATFORM = IBM_AC922 ] ; then
    TEST_COMMAND="module load $COMET_CUDA_MODULE; env OMP_NUM_THREADS=1 jsrun --nrs 2 --rs_per_host 1 --cpu_per_rs 32 -g 6 --tasks_per_rs 32 -X 1"
    TEST_COMMAND+=" -E LD_PRELOAD=${OLCF_SPECTRUM_MPI_ROOT}/lib/pami_451/libpami.so"
    CMAKE_CXX_FLAGS+=" -DHAVE_INT128"
  fi

  if [ $COMET_PLATFORM = DGX2 -o $COMET_PLATFORM = GPUSYS2 ] ; then
    CMAKE_CXX_FLAGS+=" -DHAVE_INT128"
  fi

  if [ $COMET_PLATFORM = EDISON ] ; then
    TEST_COMMAND="env OMP_NUM_THREADS=24 srun -n 64"
  fi

  #----------------------------------------------------------------------------
  #---Other.

  local CMAKE_EXTRA_OPTIONS=""

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
    CMAKE_EXTRA_OPTIONS+=" -DUSE_MPI:BOOL=ON"
  else
    CMAKE_EXTRA_OPTIONS+=" -DUSE_MPI:BOOL=OFF"
  fi

  if [ "$USE_CUDA" = ON ] ; then
    CMAKE_EXTRA_OPTIONS+=" -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
    CMAKE_EXTRA_OPTIONS+=" -DUSE_CUDA:BOOL=ON"
    CMAKE_EXTRA_OPTIONS+=" -DCUDA_HOST_COMPILER:STRING=$(dirname $(which $CC_serial))"
  else
    CMAKE_EXTRA_OPTIONS+=" -DUSE_CUDA:BOOL=OFF"
  fi

  if [ "$USE_HIP" = ON ] ; then
    CMAKE_EXTRA_OPTIONS+=" -DUSE_HIP:BOOL=ON"
  else
    CMAKE_EXTRA_OPTIONS+=" -DUSE_HIP:BOOL=OFF"
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
    -DGTEST_LD_FLAGS:STRING="$GTEST_LD_FLAGS" \
    -DTEST_COMMAND="$TEST_COMMAND" \
   \
    $CMAKE_EXTRA_OPTIONS \
   \
    $REPO_DIR

  ln -s $INSTALL_DIR install_dir
}

#    -DCUDA_HOST_COMPILER:STRING="$(dirname $(which $CC_serial))" \

#==============================================================================

main "$@"

#==============================================================================
# cruft.

#  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
#  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
#  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \

  #  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v$DEBUG_FLAG" \
  #  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \

#==============================================================================
