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
EXEC=$WD/../install/bin/genomics_metric

#---cd to Lustre to be able to aprun on titan.
pushd $MEMBERWORK/$ACCOUNT > /dev/null

#for compute_method in CPU GPU ; do
for compute_method in CPU GPU ; do

  #echo \
  #aprun -n1 $EXEC \
  #    --num_field 5000 --num_vector_local 6000 \
  #    --compute_method $compute_method --verbosity 1
  #aprun -n1 $EXEC \
  #    --num_field 5000 --num_vector_local 6000 \
  #    --compute_method $compute_method --verbosity 1

  echo \
  aprun -n1 $EXEC \
      --num_field 50 --num_vector_local 60 \
      --compute_method $compute_method --verbosity 1
  aprun -n1 $EXEC \
      --num_field 50 --num_vector_local 60 \
      --compute_method $compute_method --verbosity 1


done

popd > /dev/null

#==============================================================================
