#!/bin/bash
#==============================================================================
#
# Create modified MAGMA versions that are cloned and then patched
# with modified GEMM code.
#
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main
{
  local TAG
  for TAG in minproduct mgemm2 mgemm3 mgemm4 mgemm5 ; do
    ./clone_magma.sh $TAG
    ./patch_magma.sh $TAG
  done
}

#==============================================================================

main "$@"

#==============================================================================
