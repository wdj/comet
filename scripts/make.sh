#!/bin/bash -l
#==============================================================================
#
# Build CoMet code.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  # Location of this script.
  local SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  time make -j4 VERBOSE=1

  if [ $? != 0 ] ; then
    exit $?
  fi

  time make install
  exit $?
} # main

#==============================================================================

main "$@"

#==============================================================================
