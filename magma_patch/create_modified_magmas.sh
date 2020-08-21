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
  local VERSION="${1:-}"
  local TAG
  for TAG in minproduct mgemm2 mgemm3 mgemm4 mgemm5 ; do
    ./clone_magma.sh $TAG $VERSION
    ./patch_magma.sh $TAG $VERSION
  done
}

#==============================================================================

main "$@"

#==============================================================================
