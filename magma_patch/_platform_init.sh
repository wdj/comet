#==============================================================================
# Make initializations pertaining to platform to build/use CoMet on.
#==============================================================================

# Initial checks.

local CBE_="${COMET_BUILD_EXPERIMENTAL:-}"
if [ "$CBE_" != "" -a  \ "$CBE_" != "YES" -a "$CBE_" != "NO" ] ; then
  echo "Error in COMET_BUILD_EXPERIMENTAL setting." 1>&2
  exit 1
fi

# Set platform variables.

local COMET_PLATFORM=""
[[ -n "${CRAYOS_VERSION:-}" ]] && COMET_PLATFORM=CRAY_XK7 # OLCF Titan, Chester
[[ -n "${LSF_BINDIR:-}" ]] && COMET_PLATFORM=IBM_AC922 # OLCF Summit. Peak
[[ "$(uname -n)" = "dgx2-b" ]] && COMET_PLATFORM=DGX2 # ORNL DGX2
[[ "$(uname -n)" = "gpusys2" ]] && COMET_PLATFORM=GPUSYS2 # Turing GPU
[[ "${NERSC_HOST:-}" = "edison" ]] && COMET_PLATFORM=EDISON
[[ "${COMET_BUILD_EXPERIMENTAL:-}" = YES ]] && COMET_PLATFORM=EXPERIMENTAL
if [ "$COMET_PLATFORM" = "" ] ; then
  echo "Error: unknown platform." 1>&2
  exit 1
fi

# Modules.

if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
  true # skip
elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  module load cudatoolkit
  module load acml
  module load cmake
  (module list) 2>&1 | grep -v '^ *$'
elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
  module -q load gcc/6.4.0
  local COMET_CUDA_MODULE=cuda
  module -q load $COMET_CUDA_MODULE
  module -q load cmake
  (module list) 2>&1 | grep -v '^ *$'
elif [ $COMET_PLATFORM = DGX2 ] ; then
  true # skip
elif [ $COMET_PLATFORM = GPUSYS2 ] ; then
  true # skip
elif [ $COMET_PLATFORM = EDISON ] ; then
  module swap PrgEnv-intel PrgEnv-gnu
else
  echo "Unknown platform." 1>&2
  exit 1
fi

# Other.

local COMET_HOST
COMET_HOST="$(echo $(hostname -f) | \
              sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' \
                  -e 's/[.-].*//' -e 's/[0-9]*$//')"
local COMET_PLATFORM_STUB
[[ $COMET_PLATFORM = EXPERIMENTAL ]] && COMET_PLATFORM_STUB=experimental \
                                     || COMET_PLATFORM_STUB=$COMET_HOST

#==============================================================================