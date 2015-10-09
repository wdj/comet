#!/bin/bash -l
#==============================================================================
#
# Script to run test cases of code.
#
# Usage:
#
# qsub -I -Astf006 -lnodes=1 -lwalltime=2:0:0
# ...
# cd genomics_gpu
# cd ../build
# ../genomics_gpu/scripts/run.sh
#
#==============================================================================

ACCOUNT=stf006

WD=$PWD
#---cd to Lustre to be able to aprun on titan.
pushd $MEMBERWORK/$ACCOUNT > /dev/null

#for compute_method in CPU GPU ; do
for compute_method in GPU ; do

  echo \
  aprun -n1 $WD/../bin/genomics_metric \
      --num_field 5000 --num_vector_local 5000 \
      --compute_method $compute_method --verbosity 1
  aprun -n1 $WD/../bin/genomics_metric \
      --num_field 5000 --num_vector_local 5000 \
      --compute_method $compute_method --verbosity 1

done

popd > /dev/null

#==============================================================================
