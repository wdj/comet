#!/bin/bash
#------------------------------------------------------------------------------
# Helper script to perform build.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#------------------------------------------------------------------------------

# Exit immediately on error.
set -eu -o pipefail

HOST_=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//')
DIRNAME_STUB=$HOST_

for dir in build_*_$DIRNAME_STUB ; do
  pushd $dir
  ../genomics_gpu/scripts/make.sh 2>&1 | tee out_make.sh
  if [ $? != 0 ] ; then
    exit 1
  fi
  popd
done

#------------------------------------------------------------------------------
