#!/bin/bash
#==============================================================================
#
# Configure ALL CoMet versions prior to build.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#
# Relevant input variables:
#
# OLCF_PROJECT - OLCF project ID
# COMET_BUILD_EXPERIMENTAL -
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function script_dir
{
  local TARGET_FILE="$0"

  if [ $(uname -s) != "Darwin" ] ; then
    echo $(dirname "$(readlink -f "$TARGET_FILE")")
    return
  fi

  cd $(dirname $TARGET_FILE)
  TARGET_FILE=$(basename $TARGET_FILE)
  # Iterate down a (possible) chain of symlinks
  while [ -L "$TARGET_FILE" ] ; do
    TARGET_FILE=$(readlink $TARGET_FILE)
    cd $(dirname $TARGET_FILE)
    TARGET_FILE=$(basename $TARGET_FILE)
  done
  # Compute the canonicalized name by finding the physical path
  # for the directory we're in and appending the target file.
  local PHYS_DIR=$(pwd -P)
  local RESULT=$PHYS_DIR/$TARGET_FILE
  echo $(dirname $RESULT)
}

#==============================================================================

function configure_1case
{
  local DO_BUILD="$1"

  local BUILD_STUB=""
  [[ $FP_PRECISION = SINGLE ]] && BUILD_STUB+="single_"
  [[ $BUILD_TYPE = Debug ]] && BUILD_STUB+="test" || BUILD_STUB+="release"
  [[ $USE_MPI = OFF ]] && BUILD_STUB+="_nompi"

  local BUILD_DIR=build_${BUILD_STUB}_$COMET_PLATFORM_STUB
  local INSTALL_DIR=$INSTALLS_DIR/install_${BUILD_STUB}_$COMET_PLATFORM_STUB

  rm -rf $BUILD_DIR # Clean out any previous build files
  rm -rf $INSTALL_DIR

  if [ $DO_BUILD = ON ] ; then

    printf -- '-%.0s' {1..79}; echo ""
    echo COMET_PLATFORM_STUB=$COMET_PLATFORM_STUB INSTALLS_DIR=$INSTALLS_DIR \
      BUILD_TYPE=$BUILD_TYPE TESTING=$TESTING USE_MPI=$USE_MPI \
      FP_PRECISION=$FP_PRECISION

    echo "Creating $BUILD_DIR ..."
    mkdir $BUILD_DIR
    pushd $BUILD_DIR

    # Link to common MAGMA build if available.
    local MAGMA_BUILD_DIR=../magma_build_$COMET_PLATFORM_STUB
    if [ -e $MAGMA_BUILD_DIR ] ; then
      ln -s $MAGMA_BUILD_DIR magma_patch
    fi

    echo "Executing cmake.sh script ..."
    pwd
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=$BUILD_TYPE TESTING=$TESTING \
        USE_MPI=$USE_MPI FP_PRECISION=$FP_PRECISION \
        ../genomics_gpu/scripts/cmake.sh

    # Move magma build to location for common use for different builds.
    if [ -e magma_patch -a ! -e $MAGMA_BUILD_DIR ] ; then
      mv magma_patch $MAGMA_BUILD_DIR # share common MAGMA build
      ln -s          $MAGMA_BUILD_DIR magma_patch
    fi

    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .

    printf -- '-%.0s' {1..79}; echo ""

  fi # DO_BUILD
} # configure_1case

#==============================================================================

function main
{
  # Initial checks.
  local CBE_="${COMET_BUILD_EXPERIMENTAL:-}"
  if [ "$CBE_" != "" -a  \ "$CBE_" != "ON" -a "$CBE_" != "OFF" ] ; then
    echo "${0##*/}: Error in COMET_BUILD_EXPERIMENTAL setting." 1>&2
    exit 1
  fi

  # Location of this script.
  local SCRIPT_DIR=$(script_dir)
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  [[ -z "${OLCF_PROJECT:-}" ]] && local OLCF_PROJECT=stf006

  # Set directory for all installs.
  if [ $COMET_PLATFORM = EXPERIMENTAL ] ; then
    true # skip
  elif [ $COMET_PLATFORM = CRAY_XK7 ] ; then
    local INSTALLS_DIR=/lustre/atlas/scratch/$(whoami)/$OLCF_PROJECT/comet_work
  elif [ $COMET_PLATFORM = IBM_AC922 ] ; then
    local INSTALLS_DIR=$MEMBERWORK/$OLCF_PROJECT/comet_work
    #local INSTALLS_DIR=/gpfs/alpine/$OLCF_PROJECT/scratch/$(whoami)/comet
  elif [ $COMET_PLATFORM = EDISON ] ; then
    local INSTALLS_DIR="$SCRATCH/comet"
  else
    local INSTALLS_DIR="$PWD/installs"
    #echo "${0##*/}: Unknown platform." 1>&2
    #exit 1
  fi
  mkdir -p "$INSTALLS_DIR"

  # Prepare for builds.
  export COMET_PLATFORM_STUB=$COMET_PLATFORM_STUB INSTALLS_DIR=$INSTALLS_DIR

  #----------------------------------------------------------------------------
  # Build: test / double precision / mpi case.

  local DO_BUILD=ON # OFF
  [[ $COMET_CAN_USE_MPI = OFF ]] && DO_BUILD=OFF
  export BUILD_TYPE=Debug TESTING=ON USE_MPI=ON FP_PRECISION=DOUBLE
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build: test / double precision / nompi case.

  local DO_BUILD=ON # OFF
  export BUILD_TYPE=Debug TESTING=ON USE_MPI=OFF FP_PRECISION=DOUBLE
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build: test / single precision / mpi case.

  local DO_BUILD=ON # OFF
  [[ $COMET_CAN_USE_MPI = OFF ]] && DO_BUILD=OFF
  export BUILD_TYPE=Debug TESTING=ON USE_MPI=ON FP_PRECISION=SINGLE
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build: test / single precision / nompi case.

  local DO_BUILD=ON # OFF
  export BUILD_TYPE=Debug TESTING=ON USE_MPI=OFF FP_PRECISION=SINGLE
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build: release / double precision / mpi case.

  local DO_BUILD=ON # OFF
  [[ $COMET_CAN_USE_MPI = OFF ]] && DO_BUILD=OFF
  export BUILD_TYPE=Release TESTING=OFF USE_MPI=ON FP_PRECISION=DOUBLE
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build: release / double precision / nompi case.

  local DO_BUILD=ON # OFF
  export BUILD_TYPE=Release TESTING=OFF USE_MPI=OFF FP_PRECISION=DOUBLE
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build release / single precision / mpi case.

  local DO_BUILD=ON # OFF
  [[ $COMET_CAN_USE_MPI = OFF ]] && DO_BUILD=OFF
  export BUILD_TYPE=Release TESTING=OFF USE_MPI=ON FP_PRECISION=SINGLE 
  configure_1case $DO_BUILD

  #----------------------------------------------------------------------------
  # Build release / single precision / nompi case.

  local DO_BUILD=ON # OFF
  export BUILD_TYPE=Release TESTING=OFF USE_MPI=OFF FP_PRECISION=SINGLE 
  configure_1case $DO_BUILD
} # main

#==============================================================================

main "$@"

#==============================================================================
