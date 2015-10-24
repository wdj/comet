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
EXEC=$WD/../install_debug/bin/genomics_metric

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

#  echo \
#  aprun -n1 $EXEC \
#      --num_field 50 --num_vector_local 60 --num_way 3 \
#      --compute_method $compute_method --verbosity 1
#  aprun -n1 $EXEC \
#      --num_field 50 --num_vector_local 60 --num_way 3 \
#      --compute_method $compute_method --verbosity 1

#  for i in {1..2} ; do
#  aprun -n$i $EXEC \
#      --num_field 1 --num_vector_local $(( 6 / $i )) --num_way 2 \
#      --compute_method $compute_method --verbosity 2 --all2all yes
#  done

  #aprun -n1 -N1 $EXEC \
  #    --num_field 1 --num_vector_local 4 --all2all yes \
  #    --compute_method $compute_method --verbosity 2
  for n in 1 3 ; do
  aprun -n$n -N1 $EXEC \
      --num_field 10 --num_vector_local $(( 24 / $n )) --all2all yes \
      --compute_method $compute_method --verbosity 1
  done
  #aprun -n4 -N1 $EXEC \
  #    --num_field 1 --num_vector_local 1 --all2all yes \
  #    --compute_method $compute_method --verbosity 2

done

popd > /dev/null

#==============================================================================
