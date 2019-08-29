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

function main
{
  # Location of this script.
  local SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  . $SCRIPT_DIR/_platform_init.sh

  local DIR
  for DIR in build_*_$COMET_PLATFORM_STUB ; do
    pushd $DIR
    $SCRIPT_DIR/make.sh 2>&1 | tee out_make.sh
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
