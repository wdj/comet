#!/bin/bash

# Executable location
EXE=${COMET_SRC}/comet_work_gcc8.3/build_test_nompi_gpusys/genomics_metric

# Problem sizes
# Very Short tests
#nfields=2048
nfields=1024
#nfields=512
#nfields=256
#nfields=128

nvectors=256
#nvectors=16
#nvectors=8
#nvectors=4
#nvectors=2

# Small test - about 2s
#nfields=20480
#nvectors=2048

# Medium test - about 6s
#nfields=51200
#nvectors=5120

# Large test
#nfields=81920
#nvectors=8192

# Test options
verbose=1 # 1=minimal, 2=results, 3=all

testtype=0 # Run performance tests by default
if [ $# -gt 0 ]
then
testtype=$1
fi

# Output run options
#echo "Run options:"
#$EXE
#exit 0

if [ $testtype -eq 0 ]
then
    # Performance options
    RUNOPTS="--checksum no"
    echo "Running performance tests"
else
    # Accuracy options
    RUNOPTS="--nopreflight --print_details yes"
    echo "Running accuracy tests"
fi

PROBOPTS="--num_field $nfields --num_vector $nvectors --metric_type duo --sparse yes --compute_method GPU --all2all yes --num_way 2 --num_tc_steps 1 --verbosity $verbose"
echo "Running tests using $PROBOPTS $RUNOPTS"

# Original routine
echo -e "\n\nRunning original Magma duo CoMet GEMM"
time $EXE $PROBOPTS $RUNOPTS --num_kernel 0

#echo -e "\n\nRunning simple WMMA 1-bit duo CoMet GEMM"
#time $EXE $PROBOPTS $RUNOPTS --num_kernel 1

# Extremely slow
#echo -e "\n\nRunning very simple CoMet 1-bit xor duo general GEMM"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 0

# Fastest 1-bit WMMA GEMM tested so far
#echo -e "\n\nRunning simple tensor core 1-bit xor duo general GEMM"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 1

# Needs to be more fully optimized to actually be faster
#echo -e "\n\nRunning simple shared memory tensor core 1-bit xor duo general GEMM"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 2

#echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 256x128"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 10

# Fastest Cutlass device-level kernel on gpusys2
echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 128x256"
time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 11

#echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 128x128"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 12

#echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 128x64"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 13

#echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 64x128"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 14

#echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 64x64"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 15

#echo -e "\n\nRunning Cutlass device-level tensor core 1-bit xor duo GEMM 64x64"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 16

# Doesn't work currently, need to rearrange output array differently
#echo -e "\n\nRunning simple CoMet int xor duo GEMM"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 20

#echo -e "\n\nRunning simple CoMet xor duo GEMM"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 21

# Simple tensor core CoMet GEMM kernel
echo -e "\n\nRunning simple tensor core CoMet xor duo GEMM"
time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 22

# In progress - Optimized tensor core CoMet GEMM kernel
#echo -e "\n\nRunning optimized tensor core CoMet xor duo GEMM"
#time $EXE $PROBOPTS $RUNOPTS --tc 5 --num_kernel 23

