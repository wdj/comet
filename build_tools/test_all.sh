#!/bin/bash
#------------------------------------------------------------------------------
# Helper script to run code tests.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#------------------------------------------------------------------------------

# Exit immediately on error.
set -eu -o pipefail

# qsub -Astf006 -lnodes=4 -lwalltime=2:0:0 -I
# bsub -P stf006 -Is -nnodes 2 -alloc_flags gpumps -W 120 /bin/bash

HOST_=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//')
DIRNAME_STUB=$HOST_

if [ -n "${CRAYOS_VERSION:-}" ] ; then
  # For Titan or Chester
  if [ "$PE_ENV" = "PGI" ] ; then
    module unload PrgEnv-pgi
  fi
  module load PrgEnv-gnu
  module load cudatoolkit
else
  # For Summit or Peak
  module -q load gcc/6.4.0
  #module -q load cuda/9.1.85
  module -q load cuda
fi

dirs="build_test_$DIRNAME_STUB build_single_test_$DIRNAME_STUB"
#dirs="build_test_$DIRNAME_STUB"

for i in $dirs ; do
  echo "===================="
  echo $i
  echo "===================="
  pushd $i
  time make test ARGS=-V 2>&1 | tee out_test.txt
  if [ $? != 0 ] ; then
    exit 1
  fi
  popd
done

echo "-------------------------------------------------------------------------------"
for i in $dirs ; do
  grep -H fail $i/out_test.txt
done
echo "-------------------------------------------------------------------------------"

#------------------------------------------------------------------------------
