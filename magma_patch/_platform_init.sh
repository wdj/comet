#==============================================================================
# Build initializations pertaining to platform of CoMet build.
# This script should not be used directly but is sourced by other scripts.
#==============================================================================

#------------------------------------------------------------------------------
# Initial checks.

local CBE_="${COMET_BUILD_EXPERIMENTAL:-}"
if [ "$CBE_" != "" -a  \ "$CBE_" != "YES" -a "$CBE_" != "NO" ] ; then
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
[[ "${COMET_BUILD_EXPERIMENTAL:-}" = YES ]] && COMET_PLATFORM=EXPERIMENTAL
[[ "$COMET_HOST" = "lyra" ]] && COMET_PLATFORM=LYRA # ORNL AMD GPU system
[[ "$(uname -n)" = "hal9006" ]] && COMET_PLATFORM=AMDINTERNAL # AMD internal GPU system
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

#--------------------
if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
#--------------------

  true # skip

#--------------------
elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
#--------------------

  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  module load cudatoolkit
  module load acml
  module load cmake
  (module list) 2>&1 | grep -v '^ *$'

  local COMET_C_COMPILER=$(which cc)
  local COMET_CXX_COMPILER=$(which CC)
  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS="-march=bdver1"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_CUDA=ON

  local COMET_CUDA_COMPILE_OPTS="$CRAY_CUDATOOLKIT_INCLUDE_OPTS"
  local COMET_CUDA_LINK_OPTS="$CRAY_CUDATOOLKIT_POST_LINK_OPTS"
  COMET_CUDA_LINK_OPTS+=" -lcublas -lcudart"

  local USE_MAGMA=ON

  local USE_INT128=ON

  local COMET_TEST_COMMAND="env CRAY_CUDA_PROXY=1 OMP_NUM_THREADS=16 aprun -n64"

  # Needed for MAGMA.
  local COMET_EXTRA_LINK_OPTS="-Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib"
  COMET_EXTRA_LINK_OPTS+=" -Wl,-rpath=/opt/acml/5.3.1/gfortran64_mp/lib"

#--------------------
elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
#--------------------

  module -q load gcc/6.4.0
  module -q load cuda
  module -q load cmake
  module -q load essl
  (module list) 2>&1 | grep -v '^ *$'

  local COMET_C_COMPILER=$(which mpicc)
  local COMET_CXX_COMPILER=$(which mpiCC)
  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS="-mcpu=power9 -mtune=power9"
  COMET_EXTRA_COMPILE_OPTS+=" -mcmodel=large -m64"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_CUDA=ON

  local CUDA_ROOT="$OLCF_CUDA_ROOT"
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/targets/ppc64le-linux/lib"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/targets/ppc64le-linux/lib"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -lcublas -lcudart"

  local USE_MAGMA=ON

  local USE_INT128=ON

  local USE_CPUBLAS=ON
  local COMET_CPUBLAS_COMPILE_OPTS="-I$OLCF_ESSL_ROOT/include"
  local COMET_CPUBLAS_LINK_OPTS=""
  local XLF_DIR=$(module load xl 2>/dev/null ; echo $OLCF_XLF_ROOT)/lib
  local XLF_DIR2=$(module load xl 2>/dev/null ; echo $OLCF_XL_ROOT)/lib
  COMET_CPUBLAS_LINK_OPTS+=" -L$OLCF_ESSL_ROOT/lib64"
  COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_ESSL_ROOT/lib64 -lessl"
  COMET_CPUBLAS_LINK_OPTS+=" -L$XLF_DIR -Wl,-rpath,$XLF_DIR2 -lxlf90_r"
  COMET_CPUBLAS_LINK_OPTS+=" -lxl -lxlfmath"
  COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_GCC_ROOT/lib64"

  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 jsrun --nrs 2 --rs_per_host 1"
  COMET_TEST_COMMAND+=" --cpu_per_rs 32 -g 6 --tasks_per_rs 32 -X 1"
  #COMET_TEST_COMMAND+=" -E LD_PRELOAD=${OLCF_SPECTRUM_MPI_ROOT}/lib/pami_451/libpami.so"

#--------------------
elif [ $COMET_PLATFORM = DGX2 ] ; then
#--------------------

  local COMET_C_COMPILER=$HOME/.linuxbrew/bin/gcc-6
  local COMET_CXX_COMPILER=$HOME/.linuxbrew/bin/g++-6
  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS="-std=gnu++11"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_CUDA=ON

  local CUDA_ROOT="$HOME/cuda"
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -lcublas -lcudart"

  local USE_MAGMA=ON

  local USE_INT128=ON

#--------------------
elif [ $COMET_PLATFORM = GPUSYS2 ] ; then
#--------------------

  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"

  local COMET_C_COMPILER=$(spack location --install-dir gcc)/bin/gcc
  local COMET_CXX_COMPILER=$(spack location --install-dir gcc)/bin/g++
  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_CUDA=ON
  export CUDA_ROOT=/usr/local/cuda-10.1   # FIX export
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -lcublas -lcudart"

  local USE_MAGMA=ON

  local USE_INT128=ON

#--------------------
elif [ $COMET_PLATFORM = EDISON ] ; then
#--------------------

  module swap PrgEnv-intel PrgEnv-gnu

  local COMET_C_COMPILER=$(which cc)
  local COMET_CXX_COMPILER=$(which CC)
  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_MAGMA=OFF

  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=24 srun -n 64"

#--------------------
elif [ $COMET_PLATFORM = LYRA ] ; then
#--------------------

  module load rocm
  module load hip
  #module load rocblas
  (module list) 2>&1 | grep -v '^ *$'

  export ROCM_PATH=/opt/rocm
  export ROCBLAS_PATH=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  if [ -e $ROCBLAS_PATH ] ; then
    local BLIS_PATH=$HOME/rocBLAS/extern/blis
  else
    export ROCBLAS_PATH=/opt/rocm/rocblas
  fi

  export HIP_PATH=/opt/rocm/hip

  local COMET_C_COMPILER=$(which gcc)
  local COMET_CXX_COMPILER=$(which g++)
  local USE_GCC=OFF
  local COMET_EXTRA_COMPILE_OPTS=""

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"

  local USE_MAGMA=OFF

  if [ "${BLIS_PATH:-}" != "" ] ; then
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/zen"
    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/zen"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/zen -lblis"
  fi

  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

#--------------------
elif [ $COMET_PLATFORM = AMDINTERNAL ] ; then
#--------------------

  export ROCM_PATH=/opt/rocm
  if [ -e ~/rocBLAS/build/release/rocblas-install/rocblas ] ; then
    export ROCBLAS_PATH=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  else
    export ROCBLAS_PATH=/opt/rocm/rocblas
  fi
  export HIP_PATH=/opt/rocm/hip

  local COMET_C_COMPILER=$(which gcc)
  local COMET_CXX_COMPILER=$(which g++)
  local USE_GCC=OFF
  local COMET_EXTRA_COMPILE_OPTS=""

  local USE_OPENMP=OFF

  local USE_HIP=ON

  local USE_MAGMA=OFF

  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

#--------------------
else
#--------------------

  echo "${0##*/}: Unknown platform." 1>&2
  exit 1

fi

#==============================================================================
