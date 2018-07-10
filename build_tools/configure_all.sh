#!/bin/bash
#------------------------------------------------------------------------------
# Helper script to configure all versions prior to build.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#------------------------------------------------------------------------------

# Exit immediately on error.
set -eu -o pipefail

HOST_=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//')
DIRNAME_STUB=$HOST_

PROJECT=stf006

if [ -n "${CRAYOS_VERSION:-}" ] ; then
  # For Titan or Chester
  INSTALL_TOP_DIR=/lustre/atlas/scratch/$(whoami)/$PROJECT/comet
else
  # For Summit or Peak
  INSTALL_TOP_DIR=/gpfs/alpinetds/$PROJECT/scratch/$(whoami)/comet
fi

mkdir -p "$INSTALL_TOP_DIR"

#--------------------
# test build
#--------------------

if [ 1 = 1 ] ; then
  DIR_=build_test_$DIRNAME_STUB
  echo "Creating $DIR_ ..."
  mkdir -p $DIR_
  pushd $DIR_
  rm -rf *
  if [ -e ../magma_build_$DIRNAME_STUB ] ; then
    ln -s ../magma_build_$DIRNAME_STUB magma
  fi
  INSTALL_DIR=$INSTALL_TOP_DIR/install_test_$DIRNAME_STUB
  env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Debug TESTING=ON \
      ../genomics_gpu/scripts/cmake.sh
  if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
    mv magma ../magma_build_$DIRNAME_STUB
    ln -s    ../magma_build_$DIRNAME_STUB magma
  fi
  popd
  rm -f $(basename $INSTALL_DIR)
  ln -s $INSTALL_DIR .
fi

#--------------------
# test / single precision build
#--------------------

if [ 1 = 1 ] ; then
  DIR_=build_single_test_$DIRNAME_STUB
  echo "Creating $DIR_ ..."
  mkdir -p $DIR_
  pushd $DIR_
  rm -rf *
  if [ -e ../magma_build_$DIRNAME_STUB ] ; then
    ln -s ../magma_build_$DIRNAME_STUB magma
  fi
  INSTALL_DIR=$INSTALL_TOP_DIR/install_single_test_$DIRNAME_STUB
  env INSTALL_DIR=$INSTALL_DIR FP_PRECISION=SINGLE BUILD_TYPE=Debug TESTING=ON \
      ../genomics_gpu/scripts/cmake.sh
  if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
    mv magma ../magma_build_$DIRNAME_STUB
    ln -s    ../magma_build_$DIRNAME_STUB magma
  fi
  popd
  rm -f $(basename $INSTALL_DIR)
  ln -s $INSTALL_DIR .
fi

#--------------------
# release build
#--------------------

if [ 1 = 1 ] ; then
  DIR_=build_release_$DIRNAME_STUB
  echo "Creating $DIR_ ..."
  mkdir -p $DIR_
  pushd $DIR_
  rm -rf *
  if [ -e ../magma_build_$DIRNAME_STUB ] ; then
    ln -s ../magma_build_$DIRNAME_STUB magma
  fi
  INSTALL_DIR=$INSTALL_TOP_DIR/install_release_$DIRNAME_STUB
  env INSTALL_DIR=$INSTALL_DIR BUILD_TYPE=Release \
      ../genomics_gpu/scripts/cmake.sh
  if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
    mv magma ../magma_build_$DIRNAME_STUB
    ln -s    ../magma_build_$DIRNAME_STUB magma
  fi
  popd
  rm -f $(basename $INSTALL_DIR)
  ln -s $INSTALL_DIR .
fi

#--------------------
# release / single precision build
#--------------------

if [ 1 = 1 ] ; then
  DIR_=build_single_release_$DIRNAME_STUB
  echo "Creating $DIR_ ..."
  mkdir -p $DIR_
  pushd $DIR_
  rm -rf *
  if [ -e ../magma_build_$DIRNAME_STUB ] ; then
    ln -s ../magma_build_$DIRNAME_STUB magma
  fi
  INSTALL_DIR=$INSTALL_TOP_DIR/install_single_release_$DIRNAME_STUB
  env INSTALL_DIR=$INSTALL_DIR FP_PRECISION=SINGLE BUILD_TYPE=Release \
      ../genomics_gpu/scripts/cmake.sh
  if [ ! -e  ../magma_build_$DIRNAME_STUB ] ; then
    mv magma ../magma_build_$DIRNAME_STUB
    ln -s    ../magma_build_$DIRNAME_STUB magma
  fi
  popd
  rm -f $(basename $INSTALL_DIR)
  ln -s $INSTALL_DIR .
fi

#------------------------------------------------------------------------------
