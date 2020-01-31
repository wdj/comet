#==============================================================================
# Build initializations pertaining to platform of CoMet build.
# This script should not be used directly but is sourced by other scripts.
#==============================================================================

#------------------------------------------------------------------------------
# Initial checks.

local CBE_="${COMET_BUILD_EXPERIMENTAL:-}"
if [ "$CBE_" != "" -a  \ "$CBE_" != "ON" -a "$CBE_" != "OFF" ] ; then
  echo "${0##*/}: Error in COMET_BUILD_EXPERIMENTAL setting." 1>&2
  exit 1
fi

#------------------------------------------------------------------------------
# Set platform to build for.

local COMET_HOST
COMET_HOST="$(echo $(hostname -f) | \
              sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' \
                  -e 's/[.-].*//' -e 's/[0-9]*$//')"

local COMET_PLATFORM=""
[[ -n "${CRAYOS_VERSION:-}" ]] && COMET_PLATFORM=CRAY_XK7 # OLCF Titan, Chester
[[ -n "${LSF_BINDIR:-}" ]] && COMET_PLATFORM=IBM_AC922 # OLCF Summit. Peak
[[ "$(uname -n)" = "dgx2-b" ]] && COMET_PLATFORM=DGX2 # ORNL DGX2
[[ "$(uname -n)" = "gpusys2" ]] && COMET_PLATFORM=GPUSYS2 # Turing GPU
[[ "${NERSC_HOST:-}" = "edison" ]] && COMET_PLATFORM=EDISON
[[ "${COMET_BUILD_EXPERIMENTAL:-}" = ON ]] && COMET_PLATFORM=EXPERIMENTAL
[[ "$COMET_HOST" = "lyra" ]] && COMET_PLATFORM=LYRA # ORNL AMD GPU system
[[ "$(uname -n)" = "hal9006" ]] && COMET_PLATFORM=AMDINTERNAL # AMD internal GPU system
[[ "$COMET_HOST" = "wombat" ]] && COMET_PLATFORM=WOMBAT # ORNL HPE GPU system
[[ "$(uname -s)" = "Darwin" ]] && COMET_PLATFORM=MACOS
if [ "$COMET_PLATFORM" = "" ] ; then
  echo "${0##*/}: Unknown platform." 1>&2
  exit 1
fi

local COMET_PLATFORM_STUB
[[ $COMET_PLATFORM = EXPERIMENTAL ]] && COMET_PLATFORM_STUB=experimental \
                                     || COMET_PLATFORM_STUB=$COMET_HOST

#------------------------------------------------------------------------------
# Load needed modules and set platform-specific variables..

# Defaults.

local COMET_TEST_PROCS_MAX=64

#----------------------------------------
if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
#----------------------------------------

  true # skip

#----------------------------------------
elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
#----------------------------------------

  #---Modules etc.

  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  module load cudatoolkit
  module load acml
  module load cmake
  (module list) 2>&1 | grep -v '^ *$'

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(which cc)
  local COMET_CXX_COMPILER=$(which CC)
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS="-march=bdver1"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  local COMET_CUDA_COMPILE_OPTS="$CRAY_CUDATOOLKIT_INCLUDE_OPTS"
  local COMET_CUDA_LINK_OPTS="$CRAY_CUDATOOLKIT_POST_LINK_OPTS"
  COMET_CUDA_LINK_OPTS+=" -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=35
  local COMET_MAGMA_MAKE_INC=make.inc.titan
  # Needed for MAGMA.
  local COMET_EXTRA_LINK_OPTS="-Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib"
  COMET_EXTRA_LINK_OPTS+=" -Wl,-rpath=/opt/acml/5.3.1/gfortran64_mp/lib"

  local COMET_CAN_USE_MPI=ON
  local COMET_MPI_CMAKE_OPTS="-DMPI_C_COMPILER:STRING=$COMET_C_COMPILER"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_COMPILER:STRING=$COMET_CXX_COMPILER"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib"

  #---Testing.

  local COMET_TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64"

#----------------------------------------
elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
#----------------------------------------

  #---Modules etc.

  module -q load gcc/6.4.0
  module -q load cuda
  module -q load cmake
  module -q load essl
  (module list) 2>&1 | grep -v '^ *$'

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(which mpicc)
  local COMET_CXX_COMPILER=$(which mpiCC)
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS="-mcpu=power9 -mtune=power9"
  COMET_EXTRA_COMPILE_OPTS+=" -mcmodel=large -m64"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  local CUDA_ROOT="$OLCF_CUDA_ROOT"
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/targets/ppc64le-linux/lib"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/targets/ppc64le-linux/lib"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=70
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local USE_CPUBLAS=ON
  local COMET_CPUBLAS_COMPILE_OPTS="-I$OLCF_ESSL_ROOT/include"
  COMET_CUDA_CMAKE_OPTS+=' -DCUDA_NVCC_FLAGS="-DBLAS_H=\"essl.h\""'
  local COMET_CPUBLAS_LINK_OPTS=""
  local XLF_DIR=$(module load xl 2>/dev/null ; echo $OLCF_XLF_ROOT)/lib
  local XLF_DIR2=$(module load xl 2>/dev/null ; echo $OLCF_XL_ROOT)/lib
  COMET_CPUBLAS_LINK_OPTS+=" -L$OLCF_ESSL_ROOT/lib64"
  COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_ESSL_ROOT/lib64 -lessl"
  COMET_CPUBLAS_LINK_OPTS+=" -L$XLF_DIR -Wl,-rpath,$XLF_DIR2 -lxlf90_r"
  COMET_CPUBLAS_LINK_OPTS+=" -lxl -lxlfmath"
  COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_GCC_ROOT/lib64"

  local COMET_CAN_USE_MPI=ON

  #---Testing.

  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 jsrun --nrs 2 --rs_per_host 1"
  COMET_TEST_COMMAND+=" --cpu_per_rs 32 -g 6 --tasks_per_rs 32 -X 1"
  #COMET_TEST_COMMAND+=" -E LD_PRELOAD=${OLCF_SPECTRUM_MPI_ROOT}/lib/pami_451/libpami.so"

#----------------------------------------
elif [ $COMET_PLATFORM = DGX2 ] ; then
#----------------------------------------

  #---Compiler.

  local COMET_C_COMPILER=$HOME/.linuxbrew/bin/gcc-6
  local COMET_CXX_COMPILER=$HOME/.linuxbrew/bin/g++-6
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS="-std=gnu++11"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  local CUDA_ROOT="$HOME/cuda"
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=70
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local COMET_CAN_USE_MPI=OFF

#----------------------------------------
elif [ $COMET_PLATFORM = GPUSYS2 ] ; then
#----------------------------------------

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(spack location --install-dir gcc)/bin/gcc
  local COMET_CXX_COMPILER=$(spack location --install-dir gcc)/bin/g++
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  export CUDA_ROOT=/usr/local/cuda-10.1   # FIX export
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=75
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local COMET_CAN_USE_MPI=OFF

  #---Testing.

  local COMET_TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=1"

#----------------------------------------
elif [ $COMET_PLATFORM = EDISON ] ; then
#----------------------------------------

  #---Modules etc.

  module swap PrgEnv-intel PrgEnv-gnu

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(which cc)
  local COMET_CXX_COMPILER=$(which CC)
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  #---Libraries.

  local USE_MAGMA=OFF

  local COMET_CAN_USE_MPI=ON
  local COMET_MPI_CMAKE_OPTS="-DMPI_C_COMPILER:STRING=$COMET_C_COMPILER"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_COMPILER:STRING=$COMET_CXX_COMPILER"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH2_DIR/include"
  COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH2_DIR/lib"

  #---Testing.

  #local COMET_TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64"
  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=24 srun -n 64"

#----------------------------------------
elif [ $COMET_PLATFORM = LYRA ] ; then
#----------------------------------------

  #---Modules etc.

  module load cmake
  module load rocm
  module load hip
  #module load rocblas
  (module list) 2>&1 | grep -v '^ *$'

  export ROCM_PATH=/opt/rocm
  export HIP_PATH=/opt/rocm/hip
  # Use custom rocblas build if available.
  export ROCBLAS_PATH=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  if [ -e $ROCBLAS_PATH ] ; then
    local BLIS_PATH=$HOME/rocBLAS/extern/blis
  else
    export ROCBLAS_PATH=/opt/rocm/rocblas
  fi

  #---Compiler.

  local USE_GCC=OFF
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  #---Libraries.

  local USE_MAGMA=OFF

  if [ "${BLIS_PATH:-}" != "" ] ; then
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/zen"
    COMET_CUDA_CMAKE_OPTS=' -DCUDA_NVCC_FLAGS="-DBLAS_H=\"blis.h\""'
    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/zen"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/zen -lblis"
  fi

  local COMET_CAN_USE_MPI=OFF

#----------------------------------------
elif [ $COMET_PLATFORM = AMDINTERNAL ] ; then
#----------------------------------------

  #---Modules etc.

  export ROCM_PATH=/opt/rocm
  export HIP_PATH=/opt/rocm/hip
  # Use custom rocblas build if available.
  if [ -e ~/rocBLAS/build/release/rocblas-install/rocblas ] ; then
    export ROCBLAS_PATH=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  else
    export ROCBLAS_PATH=/opt/rocm/rocblas
  fi

  #---Compiler.

  local USE_GCC=OFF
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_HIP=ON
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  local USE_OPENMP=OFF

  #---Libraries.

  local USE_MAGMA=OFF

  local COMET_CAN_USE_MPI=OFF

#----------------------------------------
elif [ $COMET_PLATFORM = WOMBAT ] ; then
#----------------------------------------

  #---Modules etc.

  module load openmpi
  (module list) 2>&1 | grep -v '^ *$'

  #---Compiler.

  local USE_GCC=ON
  #local COMET_C_COMPILER=$(which gcc)
  #local COMET_CXX_COMPILER=$(which g++)
  #local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local COMET_C_COMPILER=$(which mpicc)
  local COMET_CXX_COMPILER=$(which mpiCC)
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"
  COMET_EXTRA_COMPILE_OPTS+=" -I$(dirname $(which mpiCC))/../include"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  CUDA_ROOT=/usr/local/cuda
  export PATH=${PATH}:$CUDA_ROOT/bin
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=70
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local COMET_CAN_USE_MPI=ON

#----------------------------------------
elif [ $COMET_PLATFORM = MACOS ] ; then
#----------------------------------------

  #---Compiler.

  local USE_GCC=OFF
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++17"

  local USE_OPENMP=ON

  #---Libraries.

  local USE_CUDA=OFF

  local USE_MAGMA=OFF

  local COMET_CAN_USE_MPI=OFF

  #---Testing.

  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"

#----------------------------------------
else
#----------------------------------------

  echo "${0##*/}: Unknown platform." 1>&2
  exit 1

fi

#==============================================================================
