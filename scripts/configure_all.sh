#!/bin/bash
#==============================================================================
#
# Configure ALL versions prior to build.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#
# Relevant input variables:
#
# OLCF_PROJECT - OLCF project ID
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  local IS_CRAY_XK7 # OLCF Titan or Chester
  [[ -n "${CRAYOS_VERSION:-}" ]] && IS_CRAY_XK7="YES" || IS_CRAY_XK7="NO"
  local IS_IBM_AC922 # OLCF Summit or Peak
  [[ -n "${LSF_BINDIR:-}" ]] && IS_IBM_AC922="YES" || IS_IBM_AC922="NO"

  if [ -z "${OLCF_PROJECT:-}" ] ; then
    local OLCF_PROJECT=stf006
  fi

  local HOST_
  HOST_=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//')
  local DIRNAME_STUB=$HOST_

  if [ $IS_CRAY_XK7 = YES ] ; then
    local INSTALLS_DIR=/lustre/atlas/scratch/$(whoami)/$OLCF_PROJECT/comet
  elif [ $IS_IBM_AC922 = YES ] ; then
    local INSTALLS_DIR=/gpfs/alpine/$OLCF_PROJECT/scratch/$(whoami)/comet
  else
    echo "Unknown platform." 1>&2
    exit 1
  fi
  mkdir -p "$INSTALLS_DIR"

  #----------------------------------------------------------------------------
  # test / double precision build

  if [ 1 = 1 ] ; then
    local BUILD_DIR=build_test_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_test_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Debug TESTING=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
      mv magma ../magma_build_$DIRNAME_STUB # share as common MAGMA build
      ln -s    ../magma_build_$DIRNAME_STUB magma
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # test / single precision build

  if [ 1 = 1 ] ; then
    local BUILD_DIR=build_single_test_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_single_test_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR FP_PRECISION=SINGLE BUILD_TYPE=Debug TESTING=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
      mv magma ../magma_build_$DIRNAME_STUB # share as common MAGMA build
      ln -s    ../magma_build_$DIRNAME_STUB magma
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # release / double precision build

  if [ 1 = 1 ] ; then
    local BUILD_DIR=build_release_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_release_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Release \
        ../genomics_gpu/scripts/cmake.sh
    if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
      mv magma ../magma_build_$DIRNAME_STUB # share as common MAGMA build
      ln -s    ../magma_build_$DIRNAME_STUB magma
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # release / double precision / nompi build

  if [ 1 = 1 ] ; then
    local BUILD_DIR=build_release_nompi_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_release_nompi_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Release NOMPI=ON \
        ../genomics_gpu/scripts/cmake.sh
    if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
      mv magma ../magma_build_$DIRNAME_STUB # share as common MAGMA build
      ln -s    ../magma_build_$DIRNAME_STUB magma
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi

  #----------------------------------------------------------------------------
  # release / single precision build

  if [ 1 = 1 ] ; then
    local BUILD_DIR=build_single_release_$DIRNAME_STUB
    echo "Creating $BUILD_DIR ..."
    mkdir -p $BUILD_DIR
    pushd $BUILD_DIR
    rm -rf *
    if [ -e ../magma_build_$DIRNAME_STUB ] ; then
      ln -s ../magma_build_$DIRNAME_STUB magma # link to common MAGMA build
    fi
    local INSTALL_DIR=$INSTALLS_DIR/install_single_release_$DIRNAME_STUB
    env INSTALL_DIR=$INSTALL_DIR FP_PRECISION=SINGLE BUILD_TYPE=Release \
        ../genomics_gpu/scripts/cmake.sh
    if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
      mv magma ../magma_build_$DIRNAME_STUB # share as common MAGMA build
      ln -s    ../magma_build_$DIRNAME_STUB magma
    fi
    popd
    rm -f $(basename $INSTALL_DIR)
    ln -s $INSTALL_DIR .
  fi
}

#==============================================================================

main "$@"

#==============================================================================
