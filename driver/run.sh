#!/bin/bash -l
#------------------------------------------------------------------------------

# qsub -I -Astf006 -lnodes=1 -lwalltime=2:0:0
# ...

# ./make.sh 2>&1 >/dev/null

WD=$PWD
#---cd to Lustre to be able to aprun on titan.
pushd $MEMBERWORK/stf006 > /dev/null
#aprun -n1 $WD/genomics_metric
echo \
aprun -n1 $WD/genomics_metric --num_field 1 --num_vector_local 2 --compute_method 0
aprun -n1 $WD/genomics_metric --num_field 1 --num_vector_local 2 --compute_method 0
popd > /dev/null

#------------------------------------------------------------------------------
