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
[[ "$COMET_HOST" = "poplar" ]] && COMET_PLATFORM=POPLAR # Cray internal system
[[ "$(uname -s)" = "Darwin" ]] && COMET_PLATFORM=MACOS
[[ "$COMET_HOST" = "va" ]] && COMET_PLATFORM=MURPHY # enclave system
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
local COMET_WERROR=ON
local COMET_USE_GTEST=ON
local COMET_USE_INT128=OFF

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

  local COMET_USE_INT128=ON

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

  #module -q load gcc/6.4.0
  module -q load gcc
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

  local COMET_USE_INT128=ON

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

  local COMET_TEST_COMMAND_PERF="env OMP_NUM_THREADS=7 jsrun --nrs 12 "
  COMET_TEST_COMMAND_PERF+="--bind packed:7 --cpu_per_rs 7 --gpu_per_rs 1 "
  COMET_TEST_COMMAND_PERF+="--rs_per_host 6 --tasks_per_rs 1 -X 1"

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

  local COMET_USE_INT128=ON

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

  if [ 0 = 1 ] ; then
    pushd ~
    git clone https://github.com/spack/spack.git
    cat <<EOF | sed -e 's/^ *//' >> ~/.bashrc
      export SPACK_ROOT=$HOME/spack
      . $SPACK_ROOT/share/spack/setup-env.sh
EOF
    . ~/.bashrc
    spack install gcc@8.3.0
    spack compiler add `spack location -i gcc@8.3.0`
    spack install cmake
    cat <<EOF | sed -e 's/^ *//' >> ~/.bashrc
      export PATH="${PATH}:/usr/local/cuda/bin"
      export PATH="${PATH}:$(spack location --install-dir cmake)/bin/"
EOF
    . ~/.bashrc
    popd
  fi

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(spack location --install-dir gcc)/bin/gcc
  local COMET_CXX_COMPILER=$(spack location --install-dir gcc)/bin/g++
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local COMET_EXTRA_COMPILE_OPTS=" -std=c++14"

  #local USE_OPENMP=ON
  local USE_OPENMP=OFF
  #local COMET_OPENMP_COMPILE_OPTS="-fopenmp"
  local COMET_OPENMP_COMPILE_OPTS=""

  #---Libraries.

  local USE_CUDA=ON
  export CUDA_ROOT=/usr/local/cuda-11.0   # FIX export
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

  #local COMET_TEST_COMMAND="env CUDA_PROXY=1"
  local COMET_TEST_COMMAND=""

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
  module load openmpi
  module load rocm
  module load hip
  #module load rocblas
  (module list) 2>&1 | grep -v '^ *$'

  #export ROCM_PATH=/opt/rocm
  #export HIP_PATH=/opt/rocm/hip
  # Use custom rocblas build if available.
  local ROCBLAS_LOCAL=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  export ROCBLAS_PATH=$ROCBLAS_LOCAL
  if [ -e $ROCBLAS_PATH ] ; then
    local BLIS_PATH=$HOME/rocBLAS/extern/blis
  else
    export ROCBLAS_PATH=/opt/rocm/rocblas
    local BLIS_PATH=$HOME/rocblas_extern/blis
  fi

  #---Compiler.

  local USE_GCC=OFF
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=hipcc
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  local USE_MAGMA=OFF

  if [ "${BLIS_PATH:-}" != "" ] ; then
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
    # ./configure --disable-threading --enable-cblas generic
  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$OLCF_OPENMPI_ROOT/include"
    local COMET_MPI_LINK_OPTS="-L$OLCF_OPENMPI_ROOT/lib -Wl,-rpath,$OLCF_OPENMPI_ROOT/lib -lmpi"
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf006
  # salloc -N1 -A stf006

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_TEST_COMMAND="module load openmpi ; env OMP_NUM_THREADS=1 mpirun --npernode 48"
  else
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

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
elif [ $COMET_PLATFORM = POPLAR ] ; then
#----------------------------------------

  local COMET_CAN_USE_MPI=OFF
  #local COMET_CAN_USE_MPI=ON

  #---Modules etc.

  module load cmake
  #module load PrgEnv-cray
  if [ $COMET_CAN_USE_MPI = ON ] ; then
    module use /home/users/twhite/share/modulefiles
    module load ompi # Trey's ompi includes rocm/3.5.0
  else
    #module load rocm-alt/2.7
    #module load rocm-alt/2.9
    #module load rocm
    #module load rocm/3.5.0
    #module load rocm-alt/3.5.0
    #module load rocm-alt/3.6.0
    module load rocm/3.6.0
  fi
  (module list) 2>&1 | grep -v '^ *$'

  export ROCM_PATH=$ROCM_PATH
  export HIP_PATH=$HIP_PATH
  # Use custom rocblas build if available.
  export ROCBLAS_PATH=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  if [ -e $ROCBLAS_PATH ] ; then
    local BLIS_PATH=$HOME/rocBLAS/extern/blis
  else
    export ROCBLAS_PATH=$ROCM_PATH/rocblas
    #local BLIS_PATH=$HOME/rocblas_extern/blis
  fi

  #---Compiler.

  local USE_GCC=OFF
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=hipcc
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  COMET_HIP_COMPILE_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"
  COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  # https://llvm.org/docs/AMDGPUUsage.html

  if [ -e $ROCM_PATH/include/gtest ] ; then
    # Poplar has gtest built-in.
    local COMET_TEST_COMPILE_OPTS=""
    local COMET_TEST_LINK_OPTS="-L$ROCM_PATH/lib64 -lgtest"
  fi

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  local USE_MAGMA=OFF

  if [ "${BLIS_PATH:-}" != "" ] ; then
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/zen"
    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/zen"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/zen -lblis"
  fi

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #local MPI_HOME=$(echo $PATH | sed 's,\(^\|.*:\)\([^:]*mvapich2[^:]*\)/bin.*,\2,')
    local COMET_MPI_COMPILE_OPTS="-I$MPI_HOME/include"
    local COMET_MPI_LINK_OPTS="-L$MPI_HOME/lib -Wl,-rpath,$MPI_HOME/lib -lmpi"
    #local COMET_MPI_CMAKE_OPTS="-DMPI_C:STRING=$COMET_C_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$MPI_HOME/include"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=\"-L$MPI_HOME/lib -lmpi\""
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX:STRING=$COMET_CXX_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$MPI_HOME/include"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=\"-L$MPI_HOME/lib -lmpi\""
  fi

  # local COMET_MPI_COMPILE_OPTS="-I$OLCF_OPENMPI_ROOT/include"
  # local COMET_MPI_LINK_OPTS="-L$OLCF_OPENMPI_ROOT/lib -Wl,-rpath,$OLCF_OPENMPI_ROOT/lib -lmpi"

  #---Testing.

  COMET_USE_GTEST=OFF

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #XXX salloc -N2 -A stf006 $SHELL
    #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
    # salloc -N1
    # salloc -N1 -pamdMI60
    # salloc -N1 -pamdMI100 --reservation=maintenance
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64 --cpu-bind=map_ldom:1 --mem-bind=local"
  else
    # salloc -N1
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1 --cpu-bind=map_ldom:1 --mem-bind=local"
  fi
  #local COMET_TEST_COMMAND="module load openmpi ; env OMP_NUM_THREADS=2 mpirun --npernode 48"
  #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -N 2 --ntasks-per-node=32"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
elif [ $COMET_PLATFORM = WOMBAT ] ; then
#----------------------------------------

  #---Modules etc.

  module load openmpi
  module load cuda
  (module list) 2>&1 | grep -v '^ *$'

  #---Compiler.

  local USE_GCC=ON
  #local COMET_C_COMPILER=$(which gcc)
  #local COMET_CXX_COMPILER=$(which g++)
  #local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local USE_PGI=NO
  if [ $USE_PGI = YES ] ; then
    local COMET_C_COMPILER=$(which pgcc)
    local COMET_CXX_COMPILER=$(which pgc++)
  else
    local COMET_C_COMPILER=$(which mpicc)
    local COMET_CXX_COMPILER=$(which mpiCC)
  fi
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"
  COMET_EXTRA_COMPILE_OPTS+=" -I$(dirname $(which mpiCC))/../include"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

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

  #local COMET_CAN_USE_MPI=ON
  local COMET_CAN_USE_MPI=OFF

  #---Testing.

  #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=224"
  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"

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
elif [ $COMET_PLATFORM = MURPHY ] ; then
#----------------------------------------

  local CONDA_PREFIX=$HOME/.conda/envs/comet_build

  if [ 0 = 1 ] ; then
  #if [ ! -e $CONDA_PREFIX/g++ ] ; then

    # NOTE: Prerequisite: need to first install conda.

    conda create --name comet_build
    conda activate comet_build
    conda install gxx_linux-64=7.3.0
    conda install openmpi
    conda deactivate

    pushd $HOME
    sh /software/source/cuda/v.10.0/cuda_10.0.130_410.48_linux.run \
      --toolkit --toolkitpath=$HOME/cuda-10.0 --silent
    popd

    pushd $HOME/cuda-10.0/bin
    ln -s $CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc gcc
    ln -s $CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-g++ g++
    popd

    pushd $HOME/.conda/envs/comet_build/bin
    ln -s x86_64-conda_cos6-linux-gnu-gcc gcc
    ln -s x86_64-conda_cos6-linux-gnu-g++ g++
    popd

  fi

  module load apps/cmake

  #---Compiler.

  local USE_GCC=ON
  local GCC_VERSION=7.3.0
  local COMET_C_COMPILER=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc
  local COMET_CXX_COMPILER=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-g++
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"
  #local COMET_EXTRA_COMPILE_OPTS+=" -I$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/sysroot/usr/include"

  #export C_INCLUDE_PATH="$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/sysroot/usr/include:$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/include/c++/${GCC_VERSION}:$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/include/c++/${GCC_VERSION}/tr1"
  export C_INCLUDE_PATH="$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/include/c++/${GCC_VERSION}:$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/sysroot/usr/include"
  export CPLUS_INCLUDE_PATH="$C_INCLUDE_PATH"
  export LIBRARY_PATH="$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib:$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/sysroot/lib"
  export PATH="$CONDA_PREFIX/bin:${PATH}"

  local USE_OPENMP=OFF # currently not available via conda.
  #local USE_OPENMP=ON
  #local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  export CUDA_ROOT=$HOME/cuda-10.0   # FIX export
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/lib64"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $COMET_CXX_SERIAL_COMPILER)
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_TOOLKIT_ROOT_DIR=$CUDA_ROOT"

  export PATH="${PATH}:$CUDA_ROOT/bin"

  local USE_MAGMA=ON
  #local USE_MAGMA=OFF
  local COMET_MAGMA_GPU_ARCH=37
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #local MPI_HOME=/usr/mpi/gcc/openmpi-4.0.3rc4 # NOTE this doesn't work
    local MPI_HOME=$CONDA_PREFIX
    local COMET_MPI_COMPILE_OPTS="-I$MPI_HOME/include"
    local COMET_MPI_LINK_OPTS="-L$MPI_HOME/lib -Wl,-rpath,$MPI_HOME/lib -lmpi_cxx -lmpi"
    local COMET_MPI_CMAKE_OPTS="-DMPI_C_COMPILER:STRING=$COMET_C_COMPILER"
    COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$MPI_HOME/include"
    COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$MPI_HOME/lib64"
    COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_COMPILER:STRING=$COMET_CXX_COMPILER"
    COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$MPI_HOME/include"
    COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$MPI_HOME/lib64"
  fi

  #---Testing.

  # 36 cores, 4 GPUS per node

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    # salloc -N 2 -n 64 -p gpu
    # salloc -N 2 -n 64 -p gpu --exclude=va-murphy-21,va-murphy-22,va-murphy-23,va-murphy-g15
    local COMET_TEST_COMMAND="env $CONDA_PREFIX/bin/mpirun -npernode 32"
  else
    # salloc -N 1 -p gpu
    local COMET_TEST_COMMAND="env srun"
  fi

#----------------------------------------
else
#----------------------------------------

  echo "${0##*/}: Unknown platform." 1>&2
  exit 1

fi

#==============================================================================
