#!/bin/bash
#==============================================================================
#
# Script to build Magma library.
#
# Usage: type the following from the Magma root directory, e.g.,
# from magma_minproduct-1.6.2:
#
# ../make.sh
#
#==============================================================================

function do_make
{
  module unload PrgEnv-pgi
  module load PrgEnv-gnu
  module load cudatoolkit
  module load acml

  time make -j8
}

#==============================================================================

do_make "$@" 2>&1 | tee out_make.txt

#==============================================================================
