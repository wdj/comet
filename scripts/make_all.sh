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
  local HOST_
  HOST_=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//')
  local DIRNAME_STUB=$HOST_

  local dir
  for dir in build_*_$DIRNAME_STUB ; do
    pushd $dir
    ../genomics_gpu/scripts/make.sh 2>&1 | tee out_make.sh
    if [ $? != 0 ] ; then
      echo "Buold failre." 1>&2
      exit 1
    fi
    popd
  done
}

#==============================================================================

main "$@"

#==============================================================================
