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

module swap PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module load acml

time make -j8 2>&1 | tee out_make.txt

#==============================================================================
