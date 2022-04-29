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

[[ "$COMET_HOST" = "node" ]] && COMET_HOST="${SLURM_SUBMIT_HOST:-}"
[[ "$COMET_HOST" = "cori" ]] && COMET_HOST="cgpu"
[[ $(hostname -f | sed -e 's/.*\.//') = "juwels" ]] && COMET_HOST="juwels"
# [[ $(echo "$COMET_HOST" | sed -e 's/.*\.//') = "jwlogin" ]] && COMET_HOST="jwlogin"
[[ "${NERSC_HOST:-}" = "perlmutter" ]] && COMET_HOST="perlmutter"
[[ "$LMOD_SYSTEM_NAME" = "frontier" ]] && COMET_HOST="frontier"

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
[[ "$COMET_HOST" = "cgpu" ]] && COMET_PLATFORM=CORI_GPU # A100s on Cori
[[ "$COMET_HOST" = "juwels" ]] && COMET_PLATFORM=JUWELS_BOOSTER # Juelich A100 system
[[ "$COMET_HOST" = "spock" ]] && COMET_PLATFORM=SPOCK # OLCF-5 EA system.
[[ "$COMET_HOST" = "birch" ]] && COMET_PLATFORM=BIRCH # OLCF-5 EA system.
[[ "$COMET_HOST" = "borg" ]] && COMET_PLATFORM=BORG # OLCF-5 MI200 EA system.
[[ "$COMET_HOST" = "bsd" ]] && COMET_PLATFORM=BSD # ORNL DGX-A100 system.
[[ "$COMET_HOST" = "perlmutter" ]] && COMET_PLATFORM=PERLMUTTER # NERSC system.
[[ "$COMET_HOST" = "bones" ]] && COMET_PLATFORM=BONES # OLCF-5 EA MI100 system.
[[ "$COMET_HOST" = "crusher" ]] && COMET_PLATFORM=CRUSHER # OLCF-5 EA MI200 system.
[[ "$COMET_HOST" = "frontier" ]] && COMET_PLATFORM=FRONTIER # OLCF-5

if [ "$COMET_PLATFORM" = "" ] ; then
  echo "${0##*/}: Unknown platform. $COMET_HOST" 1>&2
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
local COMET_COMPUTE_CAPABILITY=0

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
  local COMET_EXTRA_COMPILE_OPTS="-mcpu=power9 -mtune=power9 -Wno-maybe-uninitialized"
  COMET_EXTRA_COMPILE_OPTS+=" -mcmodel=large -m64"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  local CUDA_ROOT="$OLCF_CUDA_ROOT"
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+=" -I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+=" -I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/targets/ppc64le-linux/lib"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/targets/ppc64le-linux/lib"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"
  #COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_70,code=compute_70;-arch=sm_70"

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=70
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local USE_CPUBLAS=ON
  local COMET_CPUBLAS_COMPILE_OPTS="-I$OLCF_ESSL_ROOT/include"
  #COMET_CUDA_CMAKE_OPTS+=' -DCUDA_NVCC_FLAGS="-DBLAS_H=\"essl.h\""'
  COMET_CUDA_CMAKE_OPTS+=' -DCUDA_NVCC_FLAGS:STRING="-gencode;arch=compute_70,code=compute_70;-arch=sm_70;-DBLAS_H=\"essl.h\""'
  local COMET_CPUBLAS_LINK_OPTS=""
  local XLF_DIR=$(module load xl 2>/dev/null ; echo $OLCF_XLF_ROOT)/lib
  local XLF_DIR2=$(module load xl 2>/dev/null ; echo $OLCF_XL_ROOT)/lib
  COMET_CPUBLAS_LINK_OPTS+=" -L$OLCF_ESSL_ROOT/lib64"
  COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_ESSL_ROOT/lib64 -lessl"
  COMET_CPUBLAS_LINK_OPTS+=" -L$XLF_DIR -Wl,-rpath,$XLF_DIR2 -lxlf90_r"
  COMET_CPUBLAS_LINK_OPTS+=" -lxl -lxlfmath"
  COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_GCC_ROOT/lib64"

  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$OMPI_DIR/include"
    #local COMET_MPI_LINK_OPTS="-L$MPI_HOME/lib -Wl,-rpath,$MPI_HOME/lib -lmpi"
  fi

  #---Testing.

  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 jsrun --nrs 2 --rs_per_host 1"
  COMET_TEST_COMMAND+=" --cpu_per_rs 32 -g 6 --tasks_per_rs 32 -X 1 --smpiargs=\"-gpu\""

  local COMET_TEST_COMMAND_PERF="env OMP_NUM_THREADS=7 jsrun --nrs 12 "
  COMET_TEST_COMMAND_PERF+="--bind packed:7 --cpu_per_rs 7 --gpu_per_rs 1 "
  COMET_TEST_COMMAND_PERF+="--rs_per_host 6 --tasks_per_rs 1 -X 1 --smpiargs=\"-gpu\""

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
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_70,code=compute_70;-arch=sm_70"

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
  local COMET_C_COMPILER=$(spack location --install-dir gcc@8.3.0)/bin/gcc
  local COMET_CXX_COMPILER=$(spack location --install-dir gcc@8.3.0)/bin/g++
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local COMET_EXTRA_COMPILE_OPTS=" -std=c++14"

  local USE_OPENMP=ON
  #local USE_OPENMP=OFF
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"
  #local COMET_OPENMP_COMPILE_OPTS=""

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
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-res-usage;--ptxas-options=-v;-Xptxas;-v;-gencode;arch=compute_75,code=compute_75;-arch=sm_75"

  local USE_CUTLASS=ON
  #local COMET_CUTLASS_ARCH=Sm75
  local COMET_COMPUTE_CAPABILITY=750
  #COMET_WERROR=OFF

  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=75
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local COMET_CAN_USE_MPI=OFF

  #---Testing.

  #local COMET_TEST_COMMAND="env CUDA_PROXY=1"
  local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 "

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
  #module load hip
  #module load cray-libsci
  #module load rocblas
  (module list) 2>&1 | grep -v '^ *$'

  #export ROCM_PATH=/opt/rocm
  #export HIP_PATH=/opt/rocm/hip
  # Use custom rocblas build if available.
#  local ROCBLAS_LOCAL=$HOME/rocBLAS/build/release/rocblas-install/rocblas
#  export ROCBLAS_PATH=$ROCBLAS_LOCAL
#  if [ -e $ROCBLAS_PATH ] ; then
#    local BLIS_PATH=$HOME/rocBLAS/extern/blis
#  else
#    export ROCBLAS_PATH=/opt/rocm/rocblas
#    local BLIS_PATH=$HOME/rocblas_extern/blis
#  fi
  local ROCBLAS_PATH=$ROCM_PATH

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
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx908"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  #local USE_LAPACK=OFF
  local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
  fi

  #local USE_MAGMA=OFF
  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=gfx908

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

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
  # salloc -N2 -A stf006

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

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  #---Modules etc.

  module load cmake
  #module load PrgEnv-cray
  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #module use /home/users/twhite/share/modulefiles
    #module load ompi # Trey's ompi includes rocm/3.5.0
    #module load ompi/4.0.4-rocm-3.7
    #module load ompi/4.0.4-rocm-3.8
    #module load ompi/4.0.4-rocm-4.0
    # see https://frontier-coe.atlassian.net/wiki/spaces/FCOE/pages/109346837/Getting+Started+Guide#MPI-with-GPUs
    module load gcc/8.1.0
    module load rocm/4.1.0
    module use /home/groups/coegroup/share/coe/modulefiles
    #module load ompi/4.1.0/gnu/rocm/4.1.0
    module load ompi/4.1.0/gnu/rocm
  else
    #module load rocm-alt/2.7
    #module load rocm-alt/2.9
    #module load rocm
    #module load rocm/3.5.0
    #module load rocm-alt/3.5.0
    #module load rocm-alt/3.6.0
    #module load rocm/3.6.0
    #module load gcc/8.1.0
    #module load rocm/3.7.0
    module load gcc/8.1.0
    #module load rocm/3.8.0
    #module load rocm/4.0.0
    module load rocm/4.1.0
  fi
  (module list) 2>&1 | grep -v '^ *$'

  export ROCM_PATH=$ROCM_PATH
  export HIP_PATH=$HIP_PATH
  export ROCBLAS_PATH=$ROCM_PATH/rocblas
  ## Use custom rocblas build if available.
  #export ROCBLAS_PATH=$HOME/rocBLAS/build/release/rocblas-install/rocblas
  #if [ -e $ROCBLAS_PATH ] ; then
  #  local BLIS_PATH=$HOME/rocBLAS/extern/blis
  #else
  #  export ROCBLAS_PATH=$ROCM_PATH/rocblas
  #  #local BLIS_PATH=$HOME/rocblas_extern/blis
  #fi

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
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -Wl,-rpath,$ROCBLAS_PATH/lib -lrocblas"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -Wl,-rpath,$ROCM_PATH/lib -L$ROCM_PATH/hipsparse/lib -Wl,-rpath,$ROCM_PATH/hipsparse/lib -lrocsparse"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -Wl,-rpath,$ROCM_PATH/lib -L$ROCM_PATH/hipsparse/lib -Wl,-rpath,$ROCM_PATH/hipsparse/lib -L$ROCM_PATH/hipblas/lib -Wl,-rpath,$ROCM_PATH/hipblas/lib  $ROCM_PATH/lib/librocsparse.so  $ROCM_PATH/lib/libhipsparse.so $ROCM_PATH/hipblas/lib/libhipblas.so"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"
  COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  # https://llvm.org/docs/AMDGPUUsage.html

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx908"

  if [ -e $ROCM_PATH/include/gtest ] ; then
    # Poplar has gtest built-in.
    local COMET_TEST_COMPILE_OPTS=""
    local COMET_TEST_LINK_OPTS="-L$ROCM_PATH/lib64 -lgtest"
  fi

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  local USE_BLIS=OFF

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  fi

  #local USE_LAPACK=OFF
  local USE_LAPACK=ON

  if [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  fi

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/zen"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/zen"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/zen -lblis"
#  fi

  #local USE_MAGMA=OFF
  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=gfx908

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
    #XXX salloc -N1
    # salloc -N1 -pamdMI60
    #XXX salloc -N1 -pamdMI100 --reservation=maintenance
    # salloc -N1 -pamdMI100
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

  #module load gcc/git_master_5abe05b4
  module load cuda/11.2.1
  module load cmake

  #---Compiler.

  #local COMET_C_COMPILER=$(which gcc)
  #local COMET_CXX_COMPILER=$(which g++)
  #local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  #local USE_PGI=NO
  #if [ $USE_PGI = YES ] ; then
  #  local COMET_C_COMPILER=$(which pgcc)
  #  local COMET_CXX_COMPILER=$(which pgc++)
  #local USE_CLANG=ON
  local USE_CLANG=OFF
  if [ $USE_CLANG = ON ] ; then
    # SEGFAULTS
    local USE_GCC=OFF
    module load ARM_Compiler_For_HPC/20.3_TX2
    module load openmpi/4.0.5_armclang
    local COMET_C_COMPILER=$OMPI_CC
    local COMET_CXX_COMPILER=$OMPI_CXX
    local COMET_CXX_SERIAL_COMPILER=armclang++
    local COMET_EXTRA_COMPILE_OPTS=" -std=c++14"
  else
    local USE_GCC=ON
    module load gcc/10.2.0
    module load openmpi/4.0.5_gcc
    local COMET_C_COMPILER=$(which mpicc)
    local COMET_CXX_COMPILER=$(which mpiCC)
    local COMET_CXX_SERIAL_COMPILER=g++
    #local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++11"
    local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++14"
  fi
  COMET_EXTRA_COMPILE_OPTS+=" -I$(dirname $(which mpiCC))/../include"
  (module list) 2>&1 | grep -v '^ *$'

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  CUDA_ROOT=$CUDA_HOME
  export PATH=${PATH}:$CUDA_ROOT/bin
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  local COMET_CUDA_LINK_OPTS="-L$CUDA_ROOT/lib64"
  #COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/sbsa-linux/lib/ -lcublas -lcudart"
  COMET_CUDA_LINK_OPTS+=" -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  COMET_CUDA_LINK_OPTS+=" -rdynamic"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$(which $COMET_CXX_SERIAL_COMPILER)"
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-res-usage;--ptxas-options=-v;-Xptxas;-v;-gencode;arch=compute_70,code=compute_70;-arch=sm_70"

  #local USE_MAGMA=ON
  local USE_MAGMA=OFF
  local COMET_MAGMA_GPU_ARCH=70
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  #local COMET_CAN_USE_MPI=ON
  local COMET_CAN_USE_MPI=OFF

  #---Testing.

  #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=224"
  #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
  if [ $COMET_CAN_USE_MPI = ON ] ; then
    # salloc -N1 -pgpu
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
  else
    # salloc -N1 -pgpu
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1"
  fi

#----------------------------------------
elif [ $COMET_PLATFORM = MACOS ] ; then
#----------------------------------------

  #---Compiler.

  local USE_GCC=OFF
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=g++
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++17"

  local USE_OPENMP=OFF

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
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_37,code=compute_37;-arch=sm_37"

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
elif [ $COMET_PLATFORM = CORI_GPU ] ; then
#----------------------------------------

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  #---Modules etc.

  module load gcc # 8.3.0
  module load cuda/11.0.2
  module load cmake
  module load openmpi
  module list

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=g++
  #local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++17"
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++14"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  #local USE_CUDA=OFF
  local USE_CUDA=ON
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  #COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas_static -lcudart_static"
  COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  #local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  #COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_80,code=compute_80;-arch=sm_80"

  local USE_CUTLASS=ON
  #local COMET_CUTLASS_ARCH=Sm80
  local COMET_COMPUTE_CAPABILITY=800
  #COMET_WERROR=OFF

  #local USE_MAGMA=OFF
  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=80
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$OPENMPI_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$OPENMPI_DIR/lib -Wl,-rpath=$OPENMPI_DIR/lib -lmpi"
  fi

  #---Testing.

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    # salloc -C dgx -N 1 --ntasks-per-node=64 --cpus-per-task=1 -G 8 -t 240 -A m1759
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=16 OMP_PROC_BIND=spread OMP_PLACES=cores srun -n 64 -G 1"
  else
    # salloc -C dgx -N 1 --ntasks-per-node=1 --cpus-per-task=16 -G 8 -t 240 -A m1759
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=16 OMP_PROC_BIND=spread OMP_PLACES=cores srun -n 1"
  fi

#----------------------------------------
elif [ $COMET_PLATFORM = JUWELS_BOOSTER ] ; then
#----------------------------------------

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  #---Modules etc.

  module load GCC # 9.3.0

  # WAY 1: standard CUDA build.
  #module load CUDA # 11.0.2
  #if [ $COMET_CAN_USE_MPI = ON ] ; then
  #  module load OpenMPI # 4.1.0rc1
  #fi

  # WAY 2: custom CUDA build.
  #module use /p/scratch/share/cuda-share/modulefiles
  #module load CUDA # 11.2

  # WAY 3: custom CUDA and matching OpenMPI build.
  module use $OTHERSTAGES
  module load Stages/Devel-2020
  module load GCC CUDA/11.0 OpenMPI
  module load mpi-settings/CUDA
  module load CUDA/11.2

  module load CMake
  module list

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=g++
  #local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++17"
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++14"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  #local USE_CUDA=OFF
  local USE_CUDA=ON
  #local CUDA_ROOT=$EBROOTCUDA
  #local CUDA_ROOT=/p/project/gronor/joubert1/cuda-11.2
  # see https://stackoverflow.com/questions/19980412/how-to-let-cmake-find-cuda
  #export CUDA_BIN_PATH=$CUDA_ROOT
  #export PATH=${PATH}:CUDA_ROOT/bin
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  #COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas_static -lcudart_static"
  COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  #local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  #COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_80,code=compute_80;-arch=sm_80"

  local USE_CUTLASS=ON
  #local USE_CUTLASS=OFF
  #local COMET_CUTLASS_ARCH=Sm80
  local COMET_COMPUTE_CAPABILITY=800
  #COMET_WERROR=OFF

  #local USE_MAGMA=OFF
  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=80
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local OPENMPI_DIR=$EBROOTOPENMPI
    local COMET_MPI_COMPILE_OPTS="-I$OPENMPI_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$OPENMPI_DIR/lib -Wl,-rpath=$OPENMPI_DIR/lib -lmpi"
  fi

  #---Testing.

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    # salloc -N 2 --ntasks-per-node=32 --cpus-per-task=1 -G 8 -t 240 -A gronor -p booster
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -N 2 -n 64 -G 8"
  else
    # salloc -N 1 --ntasks-per-node=1 --cpus-per-task=24 -G 1 -t 240 -A gronor -p booster
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=24 srun -n 1 -G 1"
  fi

#----------------------------------------
elif [ $COMET_PLATFORM = SPOCK ] ; then
#----------------------------------------

  #---Modules etc.

  module load cmake/3.20.0
  module load rocm
  #module load hip
  #module load cray-libsci
  #module load rocblas
  (module list) 2>&1 | grep -v '^ *$'

  #export ROCM_PATH=/opt/rocm
  #export HIP_PATH=/opt/rocm/hip
  # Use custom rocblas build if available.
#  local ROCBLAS_LOCAL=$HOME/rocBLAS/build/release/rocblas-install/rocblas
#  export ROCBLAS_PATH=$ROCBLAS_LOCAL
#  if [ -e $ROCBLAS_PATH ] ; then
#    local BLIS_PATH=$HOME/rocBLAS/extern/blis
#  else
#    export ROCBLAS_PATH=/opt/rocm/rocblas
#    local BLIS_PATH=$HOME/rocblas_extern/blis
#  fi
  local ROCBLAS_PATH=$ROCM_PATH

  #---Compiler.

  local USE_GCC=OFF
  #local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_C_COMPILER=clang
  local COMET_CXX_COMPILER=hipcc
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  #COMET_HIP_COMPILE_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_COMPILE_OPTS+=" --offload-arch=gfx906 --offload-arch=gfx908"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  #COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_LINK_OPTS+=" --offload-arch=gfx906 --offload-arch=gfx908"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx908"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  #local USE_LAPACK=OFF
  local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
  fi

  #local USE_MAGMA=OFF
  local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=gfx908

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath,$CRAY_MPICH_DIR/lib -lmpi"

    #local COMET_MPI_CMAKE_OPTS="-DMPI_C_COMPILER:STRING=$COMET_C_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$CRAY_MPICH_DIR/include"
    ##COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib/libmpi.so"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_COMPILER:STRING=$COMET_CXX_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH_DIR/include"
    ##COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib/libmpi.so"

    local COMET_CMAKE_USE_MPI=OFF
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  # 12 compute nodes each with:
  # 1x64 core EPYC CPU
  # 256GB DDR4  Memory
  # 4xAMD MI100 GPUs with 32GiB HBM2 Memory per GPU
  # access to NCCS home and proejct areas
  # access to Alpine 

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf006 -t 360

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
  else
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1 --cpus-per-task=16 --ntasks-per-node=4 --gpu-bind=map_gpu:0,1,2,3"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
elif [ $COMET_PLATFORM = BIRCH ] ; then
#----------------------------------------

  #---Modules etc.

  #module load cmake
  module load rocm
  #module load hip
  #module load cray-libsci
  #module load rocblas
  (module list) 2>&1 | grep -v '^ *$'

  #export ROCM_PATH=/opt/rocm
  #export HIP_PATH=/opt/rocm/hip
  # Use custom rocblas build if available.
#  local ROCBLAS_LOCAL=$HOME/rocBLAS/build/release/rocblas-install/rocblas
#  export ROCBLAS_PATH=$ROCBLAS_LOCAL
#  if [ -e $ROCBLAS_PATH ] ; then
#    local BLIS_PATH=$HOME/rocBLAS/extern/blis
#  else
#    export ROCBLAS_PATH=/opt/rocm/rocblas
#    local BLIS_PATH=$HOME/rocblas_extern/blis
#  fi
  local ROCBLAS_PATH=$ROCM_PATH

  #---Compiler.

  local USE_GCC=OFF
  #local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_C_COMPILER=clang
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
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx908"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  #local USE_LAPACK=OFF
  local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
  fi

  #local USE_MAGMA=OFF
  local USE_MAGMA=ON

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON
  local COMET_MAGMA_GPU_ARCH=gfx908

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath,$CRAY_MPICH_DIR/lib -lmpi"

    #local COMET_MPI_CMAKE_OPTS="-DMPI_C_COMPILER:STRING=$COMET_C_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$CRAY_MPICH_DIR/include"
    ##COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib/libmpi.so"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_COMPILER:STRING=$COMET_CXX_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH_DIR/include"
    ##COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib/libmpi.so"

    local COMET_CMAKE_USE_MPI=OFF
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  # 12 compute nodes each with:
  # 1x64 core EPYC CPU
  # 256GB DDR4  Memory
  # 4xAMD MI100 GPUs with 32GiB HBM2 Memory per GPU
  # access to NCCS home and proejct areas
  # access to Alpine 

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf006

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
  else
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1 --cpus-per-task=16 --ntasks-per-node=4 --gpu-bind=map_gpu:0,1,2,3"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
elif [ $COMET_PLATFORM = BORG ] ; then
#----------------------------------------

  #---Modules etc.

  #module load cmake/3.20.2
  module load cmake
  module load rocm
  (module list) 2>&1 | grep -v '^ *$'

  local ROCBLAS_PATH=$ROCM_PATH

  #---Compiler.

  local USE_GCC=OFF
  #local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_C_COMPILER=clang
  local COMET_CXX_COMPILER=hipcc
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  #COMET_HIP_COMPILE_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_COMPILE_OPTS+=" --offload-arch=gfx90a"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  #COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_LINK_OPTS+=" --offload-arch=gfx90a"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx90a"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  #local USE_LAPACK=OFF
  local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
  fi

  local USE_MAGMA=OFF
  #local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=gfx90a

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath,$CRAY_MPICH_DIR/lib -lmpi"

    local COMET_CMAKE_USE_MPI=OFF
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf006 -t 360

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
  else
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1 --cpus-per-task=16 --ntasks-per-node=4 --gpu-bind=map_gpu:0,1,2,3"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
elif [ $COMET_PLATFORM = BSD ] ; then
#----------------------------------------

  #---Compiler.

  local COMET_C_COMPILER=gcc
  local COMET_CXX_COMPILER=g++
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER
  local USE_GCC=ON
  local COMET_EXTRA_COMPILE_OPTS="-std=c++14"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  local USE_CUDA=ON
  local CUDA_ROOT="/usr/local/cuda"
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  COMET_CUDA_LINK_OPTS+=" -L $CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_80,code=compute_80;-arch=sm_80"

  local USE_MAGMA=OFF
  local COMET_MAGMA_GPU_ARCH=80
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  local COMET_CAN_USE_MPI=OFF

#----------------------------------------
elif [ $COMET_PLATFORM = PERLMUTTER ] ; then
#----------------------------------------

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  #---Modules etc.

  module -q load PrgEnv-nvidia # PrgEnv-gnu
  #module -q load cpe-cuda
  module -q load cmake
  module list

  #---Compiler.

  local USE_GCC=ON
  local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_CXX_COMPILER=$(which g++) # presently unused
  local COMET_CXX_SERIAL_COMPILER=g++
  #local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++17"
  local COMET_EXTRA_COMPILE_OPTS=" -std=gnu++14"

  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local COMET_USE_INT128=ON

  #---Libraries.

  #local USE_CUDA=OFF
  local USE_CUDA=ON
  local CUDA_ROOT=$NVIDIA_PATH/cuda # $CUDA_HOME
  local COMET_CUDA_COMPILE_OPTS="-I$CUDA_ROOT/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/CUPTI/include"
  COMET_CUDA_COMPILE_OPTS+="-I$CUDA_ROOT/extras/Debugger/include"
  #COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas_static -lcudart_static"
  #COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/../math_libs/lib64 -Wl,-rpath=$CUDA_ROOT/../math_libs/lib64"
  COMET_CUDA_LINK_OPTS+=" -L$CUDA_ROOT/lib64 -Wl,-rpath=$CUDA_ROOT/lib64 -lcublas -lcudart"
  local COMET_CUDA_CMAKE_OPTS="-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON"
  #local _COMPILER_DIR_TMP_=$(dirname $(which $COMET_CXX_SERIAL_COMPILER))
  #COMET_CUDA_CMAKE_OPTS+=" -DCUDA_HOST_COMPILER:STRING=$_COMPILER_DIR_TMP_"
  COMET_CUDA_CMAKE_OPTS+=" -DCUDA_NVCC_FLAGS:STRING=-gencode;arch=compute_80,code=compute_80;-arch=sm_80"

  local USE_CUTLASS=ON
  #local USE_CUTLASS=OFF
  #local COMET_CUTLASS_ARCH=Sm80
  local COMET_COMPUTE_CAPABILITY=800
  #COMET_WERROR=OFF

  local USE_MAGMA=OFF
  #local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=80
  local COMET_MAGMA_MAKE_INC=make.inc.summit

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #local COMET_MPI_CMAKE_OPTS="-DMPI_C_COMPILER:STRING=$COMET_C_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_INCLUDE_PATH:STRING=$CRAY_MPICH_DIR/include"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_C_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_COMPILER:STRING=$COMET_CXX_COMPILER"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_INCLUDE_PATH:STRING=$CRAY_MPICH_DIR/include"
    #COMET_MPI_CMAKE_OPTS+=" -DMPI_CXX_LIBRARIES:STRING=$CRAY_MPICH_DIR/lib"
    #local OPENMPI_DIR=$EBROOTOPENMPI
    #local COMET_MPI_COMPILE_OPTS="-I$OPENMPI_DIR/include"
    #local COMET_MPI_LINK_OPTS="-L$OPENMPI_DIR/lib -Wl,-rpath=$OPENMPI_DIR/lib -lmpi"
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath=$CRAY_MPICH_DIR/lib -lmpi"
  fi

  #---Testing.

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    # salloc --nodes 2 --qos interactive --time 04:00:00 --constraint gpu --gpus-per-node 4 --account=m1759_g
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=4 OMP_PROC_BIND=spread OMP_PLACES=threads srun -N 2 -n 64 --gpus-per-node 4 --cpus-per-task=2"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 -n 64 --gpus-per-node 4"
  else
    # salloc -N 1 --ntasks-per-node=1 --cpus-per-task=24 -G 1 -t 240 -A gronor -p booster
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=64 srun -n 1 -G 1"
  fi

#----------------------------------------
elif [ $COMET_PLATFORM = BONES ] ; then
#----------------------------------------

  #---Modules etc.

  #module load cmake/3.20.2
  module load cmake
  module load rocm
  (module list) 2>&1 | grep -v '^ *$'

  local ROCBLAS_PATH=$ROCM_PATH

  #---Compiler.

  local USE_GCC=OFF
  #local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_C_COMPILER=clang
  local COMET_CXX_COMPILER=hipcc
  local COMET_CXX_SERIAL_COMPILER=hipcc

  local USE_OPENMP=OFF

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  #COMET_HIP_COMPILE_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_COMPILE_OPTS+=" --offload-arch=gfx908"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  #COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_LINK_OPTS+=" --offload-arch=gfx908"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx908"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  #local USE_LAPACK=OFF
  local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    local USE_CPUBLAS=ON
    local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
  fi

  local USE_MAGMA=OFF
  #local USE_MAGMA=ON
  local COMET_MAGMA_GPU_ARCH=gfx908

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath,$CRAY_MPICH_DIR/lib -lmpi"

    local COMET_CMAKE_USE_MPI=OFF
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf006 -t 360

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
  else
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1 --cpus-per-task=16 --ntasks-per-node=4 --gpu-bind=map_gpu:0,1,2,3"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
elif [ $COMET_PLATFORM = CRUSHER ] ; then
#----------------------------------------

  #---Modules etc.

  #module load cmake/3.20.2
  module load cmake
  #module load rocm
  module load rocm/4.5.2
  (module list) 2>&1 | grep -v '^ *$'

  local ROCBLAS_PATH=$ROCM_PATH

  #---Compiler.

  local USE_GCC=OFF
  #local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_C_COMPILER=clang
  local COMET_CXX_COMPILER=hipcc # amdclang
  local COMET_CXX_SERIAL_COMPILER=$COMET_CXX_COMPILER

  #local USE_OPENMP=OFF
  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  #COMET_HIP_COMPILE_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_COMPILE_OPTS+=" --offload-arch=gfx90a"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  #COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_LINK_OPTS+=" --offload-arch=gfx90a"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  # need this for amdclang
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_AMD__"

  if [ $USE_OPENMP = ON ] ; then
    #COMET_HIP_LINK_OPTS+=" $ROCM_PATH/llvm/lib-debug/libomp.so"
    COMET_HIP_LINK_OPTS+=" $ROCM_PATH/llvm/lib/libomp.so"
    COMET_HIP_LINK_OPTS+=" -Wl,-rpath,$ROCM_PATH/llvm/lib"
  fi

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx90a"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  local USE_LAPACK=OFF
  #local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    #local USE_CPUBLAS=ON
    #local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    #local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    #COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
    local USE_CPUBLAS=ON
    module load openblas
    local COMET_CPUBLAS_COMPILE_OPTS="-I$OLCF_OPENBLAS_ROOT/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$OLCF_OPENBLAS_ROOT/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_OPENBLAS_ROOT/lib -lopenblas"
  fi

  local USE_MAGMA=ON
  #local USE_MAGMA=OFF
  local COMET_MAGMA_GPU_ARCH=gfx90a
  #local COMET_MAGMA_MAKE_INC=make.inc.frontier

  if [ $USE_MAGMA = ON ] ; then
    COMET_HIP_LINK_OPTS+=" -Wl,-rpath,$ROCM_PATH/lib -lhipblas -lhipsparse -lrocsparse"
  fi

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath,$CRAY_MPICH_DIR/lib -lmpi"

    local COMET_CMAKE_USE_MPI=OFF
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf006 -t 360

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N2 -n64 --cpus-per-task=2 --ntasks-per-node=32 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N2 -n64 --cpus-per-task=2 --ntasks-per-node=32 --gpu-bind=closest --gpus-per-node=8 -u"
  else
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n1 --cpus-per-task=16 --ntasks-per-node=4 --gpu-bind=map_gpu:0,1,2,3"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -n1 --cpus-per-task=16 --ntasks-per-node=1 --gpu-bind=closest --gpus-per-node=8 -u"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
elif [ $COMET_PLATFORM = FRONTIER ] ; then
#----------------------------------------

  #---Modules etc.

  #module load cmake/3.20.2
  #module load cmake
  #module load rocm
  module load rocm/4.5.2
  (module list) 2>&1 | grep -v '^ *$'

  local ROCBLAS_PATH=$ROCM_PATH

  #---Compiler.

  local USE_GCC=OFF
  #local COMET_C_COMPILER=$(which gcc) # presently unused
  local COMET_C_COMPILER=clang
  local COMET_CXX_COMPILER=hipcc
  local COMET_CXX_SERIAL_COMPILER=hipcc

  #local USE_OPENMP=OFF
  local USE_OPENMP=ON
  local COMET_OPENMP_COMPILE_OPTS="-fopenmp"

  local USE_HIP=ON
  local COMET_HIP_COMPILE_OPTS="-I$ROCBLAS_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$ROCM_PATH/include"
  COMET_HIP_COMPILE_OPTS+=" -I$HIP_PATH/include/hip"
  COMET_HIP_COMPILE_OPTS+=" -fno-gpu-rdc -Wno-unused-command-line-argument"
  #COMET_HIP_COMPILE_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_COMPILE_OPTS+=" --offload-arch=gfx90a"
  COMET_HIP_COMPILE_OPTS+=" -Wno-c99-designator"
  COMET_HIP_COMPILE_OPTS+=" -Wno-duplicate-decl-specifier -Wno-unused-variable" # FIX this later after compiler headers fixed
  #COMET_HIP_COMPILE_OPTS+=" -DCUBLAS_V2_H_ -DHAVE_HIP"
  COMET_HIP_COMPILE_OPTS+=" -DHAVE_HIP"
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_HCC__"
  local COMET_HIP_LINK_OPTS="-L$ROCBLAS_PATH/lib -lrocblas"
  COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lrocsparse"
  #COMET_HIP_LINK_OPTS+=" --amdgpu-target=gfx906,gfx908"
  COMET_HIP_LINK_OPTS+=" --offload-arch=gfx90a"
  #COMET_HIP_LINK_OPTS+=" -L$ROCM_PATH/lib -lhip_hcc"

  # need this for amdclang
  #COMET_HIP_COMPILE_OPTS+=" -D__HIP_PLATFORM_AMD__"

  if [ $USE_OPENMP = ON ] ; then
    #COMET_HIP_LINK_OPTS+=" $ROCM_PATH/llvm/lib-debug/libomp.so"
    COMET_HIP_LINK_OPTS+=" $ROCM_PATH/llvm/lib/libomp.so"
    COMET_HIP_LINK_OPTS+=" -Wl,-rpath,$ROCM_PATH/llvm/lib"
  fi

  local COMET_HIP_CMAKE_OPTS="-DCOMET_HIP_ARCHITECTURES=gfx90a"

  COMET_WERROR=OFF

# If you have device code that calls other device code that exists only in the same translation unit then you can compile with the '-fno-gpu-rdc' option.  This forces the AMD compiler to emit device code at compile time rather than link time.  Link times can be much shorter.  Compile times can increase slightly you're probably already doing a parallel compile via `make -j`.

  #---Libraries.

  #local USE_BLIS=ON
  local USE_BLIS=OFF

  local USE_LAPACK=OFF
  #local USE_LAPACK=ON

  if [ "${USE_BLIS:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  elif [ "${USE_LAPACK:-OFF}" != OFF ] ; then
    local USE_CPUBLAS=ON
  else
    #local USE_CPUBLAS=ON
    #local COMET_CPUBLAS_COMPILE_OPTS="-I$CRAY_LIBSCI_PREFIX/include"
    #local COMET_CPUBLAS_LINK_OPTS="-L$CRAY_LIBSCI_PREFIX/lib"
    #COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$CRAY_LIBSCI_PREFIX/lib -lsci_cray"
    local USE_CPUBLAS=ON
    module load openblas
    local COMET_CPUBLAS_COMPILE_OPTS="-I$OLCF_OPENBLAS_ROOT/include"
    local COMET_CPUBLAS_LINK_OPTS="-L$OLCF_OPENBLAS_ROOT/lib"
    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$OLCF_OPENBLAS_ROOT/lib -lopenblas"
  fi

  local USE_MAGMA=ON
  #local USE_MAGMA=OFF
  local COMET_MAGMA_GPU_ARCH=gfx90a
  #local COMET_MAGMA_MAKE_INC=make.inc.frontier

  if [ $USE_MAGMA = ON ] ; then
    COMET_HIP_LINK_OPTS+=" -Wl,-rpath,$ROCM_PATH/lib -lhipblas -lhipsparse -lrocsparse"
  fi

#  if [ "${BLIS_PATH:-}" != "" ] ; then
#    local USE_CPUBLAS=ON
#    local COMET_CPUBLAS_COMPILE_OPTS="-I$BLIS_PATH/include/generic"
#    #COMET_CPUBLAS_COMPILE_OPTS+=' -include "blis.h"'
#    local COMET_CPUBLAS_LINK_OPTS="-L$BLIS_PATH/lib/generic"
#    COMET_CPUBLAS_LINK_OPTS+=" -Wl,-rpath,$BLIS_PATH/lib/generic -lblis"
#    # ./configure --disable-threading --enable-cblas generic
#  fi

  #local COMET_CAN_USE_MPI=OFF
  local COMET_CAN_USE_MPI=ON

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    local COMET_MPI_COMPILE_OPTS="-I$CRAY_MPICH_DIR/include"
    local COMET_MPI_LINK_OPTS="-L$CRAY_MPICH_DIR/lib -Wl,-rpath,$CRAY_MPICH_DIR/lib -lmpi"

    local COMET_CMAKE_USE_MPI=OFF
  fi

  #---Testing.

  #COMET_USE_GTEST=OFF

  #XXX salloc -N2 -A stf006 $SHELL
  #XXX srun -N 1 --ntasks-per-node=1 -A stf006  --pty bash
  #XXX salloc -N2 -A stf006 $SHELL
  #XXX salloc -N2 -A stf006 $SHELL
  # salloc -N1 -A stf016_frontier -t 360

  if [ $COMET_CAN_USE_MPI = ON ] ; then
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1 srun -n64"
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N2 -n64 --cpus-per-task=2 --ntasks-per-node=32 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N2 -n64 --cpus-per-task=2 --ntasks-per-node=32 --gpu-bind=closest --gpus-per-node=8 -u"
  else
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=1"
    #local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -n1 --cpus-per-task=16 --ntasks-per-node=1 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7"
    local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -n1 --cpus-per-task=16 --ntasks-per-node=1 --gpu-bind=closest --gpus-per-node=8 -u"
  fi
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 2 --ntasks-per-node=48"
  #XXX local COMET_TEST_COMMAND="env OMP_NUM_THREADS=2 srun -N 1 --ntasks-per-node=1"

#----------------------------------------
else
#----------------------------------------

  echo "${0##*/}: Unknown platform. $COMET_HOST $COMET_PLATFORM" 1>&2
  exit 1

fi

#==============================================================================
