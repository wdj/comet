#!/bin/bash
#==============================================================================
#
# Script to build Magma library.
#
# Usage: type the following from the Magma root directory, e.g.,
# magma_minproduct-1.6.2:
#
# ../make.sh
#
#==============================================================================

module swap PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module load acml

make -j8

#==============================================================================
