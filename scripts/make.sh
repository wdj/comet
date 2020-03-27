#!/bin/bash -l
#==============================================================================
#
# Build CoMet code.
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

function main
{
  # Location of this script.
  #local SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  local SCRIPT_DIR=$(script_dir)
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  time make -j4 VERBOSE=1

  if [ $? != 0 ] ; then
    exit $?
  fi

  #FIX time make install
  exit $?
} # main

#==============================================================================

main "$@"

#==============================================================================
