#!/bin/bash
#==============================================================================
#
# Build ALL versions.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
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

function main
{
  # Location of this script.
  local SCRIPT_DIR=$(script_dir)
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  if [ -z "${1:-}" ] ; then
    local DIRS="$(ls -d build_*_$COMET_PLATFORM_STUB)"
    #local DIRS="$(ls -d build_single_release_$COMET_PLATFORM_STUB)"
  else
    # Build only versions needed for debugging work.
    if [ $COMET_CAN_USE_MPI = ON ] ; then
      local DIRS="build_test_$COMET_PLATFORM_STUB \
                  build_single_test_$COMET_PLATFORM_STUB"
    else
      local DIRS="build_test_nompi_$COMET_PLATFORM_STUB \
                  build_single_test_nompi_$COMET_PLATFORM_STUB"
    fi
  fi

  local DIR
  for DIR in $DIRS ; do
    pushd $DIR
    $SCRIPT_DIR/make.sh 2>&1 | tee out_make.txt
    local MYSTATUS=$?
    if [ $MYSTATUS != 0 ] ; then
      echo "Build failure." 1>&2
      exit $MYSTATUS
    fi
    popd
  done
}

#==============================================================================

main "$@"

#==============================================================================
