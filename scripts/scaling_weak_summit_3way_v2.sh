#!/bin/bash
#==============================================================================
#
# Script for executing 3-way weak scaling study on Summit.
#
# Usage:
#
# env num_node_solve=1    bsub -P stf006 -nnodes   21 -W 60 -J CoMet-benchmark_1_1      scaling_weak_summit_3way_v2.sh
# env num_node_solve=1    bsub -P stf006 -nnodes   21 -W 60 -J CoMet-benchmark_2_4      scaling_weak_summit_3way_v2.sh
# env num_node_solve=2    bsub -P stf006 -nnodes   22 -W 60 -J CoMet-benchmark_3_12     scaling_weak_summit_3way_v2.sh
# env num_node_solve=4    bsub -P stf006 -nnodes   24 -W 60 -J CoMet-benchmark_4_20     scaling_weak_summit_3way_v2.sh
# env num_node_solve=6    bsub -P stf006 -nnodes   26 -W 60 -J CoMet-benchmark_5_35     scaling_weak_summit_3way_v2.sh
# env num_node_solve=10   bsub -P stf006 -nnodes   30 -W 60 -J CoMet-benchmark_6_60     scaling_weak_summit_3way_v2.sh
# env num_node_solve=14   bsub -P stf006 -nnodes   34 -W 60 -J CoMet-benchmark_7_84     scaling_weak_summit_3way_v2.sh
# env num_node_solve=20   bsub -P stf006 -nnodes   40 -W 60 -J CoMet-benchmark_8_120    scaling_weak_summit_3way_v2.sh
# env num_node_solve=29   bsub -P stf006 -nnodes   49 -W 60 -J CoMet-benchmark_9_171    scaling_weak_summit_3way_v2.sh
# env num_node_solve=37   bsub -P stf006 -nnodes   57 -W 60 -J CoMet-benchmark_10_220   scaling_weak_summit_3way_v2.sh
# env num_node_solve=48   bsub -P stf006 -nnodes   68 -W 60 -J CoMet-benchmark_11_286   scaling_weak_summit_3way_v2.sh
# env num_node_solve=62   bsub -P stf006 -nnodes   82 -W 60 -J CoMet-benchmark_12_372   scaling_weak_summit_3way_v2.sh
# env num_node_solve=76   bsub -P stf006 -nnodes   96 -W 60 -J CoMet-benchmark_13_455   scaling_weak_summit_3way_v2.sh
# env num_node_solve=94   bsub -P stf006 -nnodes  114 -W 60 -J CoMet-benchmark_14_560   scaling_weak_summit_3way_v2.sh
# env num_node_solve=115  bsub -P stf006 -nnodes  135 -W 60 -J CoMet-benchmark_15_690   scaling_weak_summit_3way_v2.sh
# env num_node_solve=136  bsub -P stf006 -nnodes  156 -W 60 -J CoMet-benchmark_16_816   scaling_weak_summit_3way_v2.sh
# env num_node_solve=162  bsub -P stf006 -nnodes  182 -W 60 -J CoMet-benchmark_17_969   scaling_weak_summit_3way_v2.sh
# env num_node_solve=192  bsub -P stf006 -nnodes  212 -W 60 -J CoMet-benchmark_18_1152  scaling_weak_summit_3way_v2.sh
# env num_node_solve=222  bsub -P stf006 -nnodes  242 -W 60 -J CoMet-benchmark_19_1330  scaling_weak_summit_3way_v2.sh
# env num_node_solve=257  bsub -P stf006 -nnodes  277 -W 60 -J CoMet-benchmark_20_1540  scaling_weak_summit_3way_v2.sh
# env num_node_solve=298  bsub -P stf006 -nnodes  318 -W 60 -J CoMet-benchmark_21_1785  scaling_weak_summit_3way_v2.sh
# env num_node_solve=338  bsub -P stf006 -nnodes  358 -W 60 -J CoMet-benchmark_22_2024  scaling_weak_summit_3way_v2.sh
# env num_node_solve=384  bsub -P stf006 -nnodes  404 -W 60 -J CoMet-benchmark_23_2300  scaling_weak_summit_3way_v2.sh
# env num_node_solve=436  bsub -P stf006 -nnodes  456 -W 60 -J CoMet-benchmark_24_2616  scaling_weak_summit_3way_v2.sh
# env num_node_solve=488  bsub -P stf006 -nnodes  508 -W 60 -J CoMet-benchmark_25_2925  scaling_weak_summit_3way_v2.sh
# env num_node_solve=546  bsub -P stf006 -nnodes  566 -W 60 -J CoMet-benchmark_26_3276  scaling_weak_summit_3way_v2.sh
# env num_node_solve=612  bsub -P stf006 -nnodes  632 -W 60 -J CoMet-benchmark_27_3672  scaling_weak_summit_3way_v2.sh
# env num_node_solve=677  bsub -P stf006 -nnodes  697 -W 60 -J CoMet-benchmark_28_4060  scaling_weak_summit_3way_v2.sh
# env num_node_solve=750  bsub -P stf006 -nnodes  770 -W 60 -J CoMet-benchmark_29_4495  scaling_weak_summit_3way_v2.sh
# env num_node_solve=830  bsub -P stf006 -nnodes  850 -W 60 -J CoMet-benchmark_30_4980  scaling_weak_summit_3way_v2.sh
# env num_node_solve=910  bsub -P stf006 -nnodes  930 -W 60 -J CoMet-benchmark_31_5456  scaling_weak_summit_3way_v2.sh
# env num_node_solve=998  bsub -P stf006 -nnodes 1018 -W 60 -J CoMet-benchmark_32_5984  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1095 bsub -P stf006 -nnodes 1115 -W 60 -J CoMet-benchmark_33_6567  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1190 bsub -P stf006 -nnodes 1210 -W 60 -J CoMet-benchmark_34_7140  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1295 bsub -P stf006 -nnodes 1315 -W 60 -J CoMet-benchmark_35_7770  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1410 bsub -P stf006 -nnodes 1430 -W 60 -J CoMet-benchmark_36_8460  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1524 bsub -P stf006 -nnodes 1544 -W 60 -J CoMet-benchmark_37_9139  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1647 bsub -P stf006 -nnodes 1667 -W 60 -J CoMet-benchmark_38_9880  scaling_weak_summit_3way_v2.sh
# env num_node_solve=1781 bsub -P stf006 -nnodes 1801 -W 60 -J CoMet-benchmark_39_10686 scaling_weak_summit_3way_v2.sh
# env num_node_solve=1914 bsub -P stf006 -nnodes 1934 -W 60 -J CoMet-benchmark_40_11480 scaling_weak_summit_3way_v2.sh
# env num_node_solve=2057 bsub -P stf006 -nnodes 2077 -W 60 -J CoMet-benchmark_41_12341 scaling_weak_summit_3way_v2.sh
# env num_node_solve=2212 bsub -P stf006 -nnodes 2232 -W 60 -J CoMet-benchmark_42_13272 scaling_weak_summit_3way_v2.sh
# env num_node_solve=2365 bsub -P stf006 -nnodes 2385 -W 60 -J CoMet-benchmark_43_14190 scaling_weak_summit_3way_v2.sh
# env num_node_solve=2530 bsub -P stf006 -nnodes 2550 -W 60 -J CoMet-benchmark_44_15180 scaling_weak_summit_3way_v2.sh
# env num_node_solve=2708 bsub -P stf006 -nnodes 2728 -W 60 -J CoMet-benchmark_45_16245 scaling_weak_summit_3way_v2.sh
# env num_node_solve=2883 bsub -P stf006 -nnodes 2903 -W 60 -J CoMet-benchmark_46_17296 scaling_weak_summit_3way_v2.sh
# env num_node_solve=3071 bsub -P stf006 -nnodes 3091 -W 60 -J CoMet-benchmark_47_18424 scaling_weak_summit_3way_v2.sh
# env num_node_solve=3272 bsub -P stf006 -nnodes 3292 -W 60 -J CoMet-benchmark_48_19632 scaling_weak_summit_3way_v2.sh
# env num_node_solve=3471 bsub -P stf006 -nnodes 3491 -W 60 -J CoMet-benchmark_49_20825 scaling_weak_summit_3way_v2.sh
# env num_node_solve=3684 bsub -P stf006 -nnodes 3704 -W 60 -J CoMet-benchmark_50_22100 scaling_weak_summit_3way_v2.sh
# env num_node_solve=3910 bsub -P stf006 -nnodes 3930 -W 60 -J CoMet-benchmark_51_23460 scaling_weak_summit_3way_v2.sh
# env num_node_solve=4134 bsub -P stf006 -nnodes 4154 -W 60 -J CoMet-benchmark_52_24804 scaling_weak_summit_3way_v2.sh
# env num_node_solve=4373 bsub -P stf006 -nnodes 4393 -W 60 -J CoMet-benchmark_53_26235 scaling_weak_summit_3way_v2.sh
#
# Options:
#
# export metric_type=ccc sparse=yes # default
# export metric_type=ccc sparse=yes tc=0
# export metric_type=ccc sparse=no
# export metric_type=ccc sparse=no tc=0
# export metric_type=duo sparse=yes
# export metric_type=czekanowski single=1
# export metric_type=czekanowski single=0
#
# original file: /ccs/home/joubert/proj/genomics/results/gbrev_max
# [SEE /ccs/home/joubert/genomics/results/gbrev_3way/run.sh]
#
#==============================================================================

#------------------------------------------------------------------------------
# Node counts

num_node_job=$(( ( $(echo $LSB_MCPU_HOSTS | wc -w) - 2 ) / 2 ))

[[ -z "${num_node_solve:-}" ]] && num_node_solve=$num_node_job
[[ -z "${num_node_launch:-}" ]] && num_node_launch=$num_node_job

ranks_per_node=6
load=6

# Find largest problem that will fit.

num_proc_field=1

for num_proc_vector_new in {1..100} ; do
  num_proc_repl_new=$(( ( ( $num_proc_vector_new + 1 ) * \
    ( $num_proc_vector_new + 2 ) + $load - 1 ) / $load ))
  num_proc_new=$(( $num_proc_vector_new * $num_proc_field * \
    $num_proc_repl_new ))
  if [ $num_proc_new -gt $(( $num_node_solve * $ranks_per_node )) ] ; then
    break
  fi
  num_proc_vector=$num_proc_vector_new
  num_proc_repl=$num_proc_repl_new
done

num_proc=$(( $num_proc_vector * $num_proc_field * $num_proc_repl ))

#------------------------------------------------------------------------------
# Algorithm settings

[[ "${ccc:-}" = 1 ]] && metric_type=ccc   # legacy settings
[[ "${ccc:-}" = 0 ]] && metric_type=czekanowski   # legacy settings
[[ -z "${metric_type:-}" ]] && metric_type=ccc
[[ -z "${single:-}" ]] && single=1
[[ -z "${sparse:-}" ]] && sparse=yes
[[ -z "${tc:-}" && $metric_type != czekanowski ]] && tc=1
[[ -z "${tc:-}" && $metric_type = czekanowski ]] && tc=0
[[ -z "${cksum:-}" ]] && cksum=yes # alt. cksum=no
[[ -z "${problem_type:-}" ]] && problem_type=random # alt. problem_type=analytic
[[ -z "${debug:-}" ]] && debug=0
[[ -z "${num_tc_steps:-}" ]] && num_tc_steps=4
#[[ -z "${num_tc_steps:-}" ]] && num_tc_steps=2
#[[ -z "${num_tc_steps:-}" ]] && num_tc_steps=5

#------------------------------------------------------------------------------
# Problem sizes

if [ "$metric_type" != czekanowski ] ; then

  #num_field_local=$(( 100000 * $num_tc_steps ))
  #num_field_local=$(( 100000 * $num_tc_steps ))
  num_field_local=$(( 100 * 1024 * $num_tc_steps ))

  #num_field_local=$(( 98304 * $num_tc_steps ))
  #num_field_local=$(( 2 * 98304 * $num_tc_steps ))
  #num_field_local=$(( 200000 * $num_tc_steps ))
  #num_field_local=$(( 2 * 100352 * $num_tc_steps ))
  #num_field_local=$(( 5 * 16384 * $num_tc_steps ))
  #num_field_local=$(( 98304 * $num_tc_steps ))

  #num_vector_local=$(( 992 * 6 ))
  num_vector_local=$(( 36 * 32 * 6 ))

  #num_stage=$(( ( $num_vector_local / 6 ) / 4 ))
  num_stage=$(( ( $num_vector_local / 6 ) / 16 ))

elif [ "$single" = 1 ] ; then
  num_field_local=$(( 7500 * $num_tc_steps ))
  num_vector_local=$(( 768 * 6 ))
  num_stage=$(( 18 ))
else
  num_field_local=$(( 7500 * $num_tc_steps ))
  num_vector_local=$(( 480 * 6 ))
  num_stage=$(( 6 ))
fi

num_vector=$(( $num_vector_local * $num_proc_vector ))

# Compute one phase of results out of possibly many phases
#[[ -z "$num_phase_ratio" ]] && num_phase_ratio=$(( 30 * $num_proc_repl ))
#[[ -z "$num_phase" ]] && num_phase=$(( ( $num_proc_vector + $num_phase_ratio - 1 ) / $num_phase_ratio ))
num_phase=1
[[ -z "$phase_min" ]] && phase_min=$(( $num_phase - 4 ))
[[ $phase_min < 0 ]] && phase_min=0
[[ -z "$phase_max" ]] && phase_max=$(( $num_phase - 1 ))

#------------------------------------------------------------------------------
# Execution settings

module -q load gcc
module -q load cuda

host=$(echo $(hostname -f) | sed -e 's/^login[0-9]\.//' -e 's/^batch[0-9]\.//' -e 's/[.-].*//' -e 's/[0-9]*$//')
OLCF_PROJECT="$(echo $LSB_PROJECT_NAME | tr A-Z a-z)"
#INSTALLS_DIR=/gpfs/alpine/$OLCF_PROJECT/scratch/$(whoami)/comet
INSTALLS_DIR=$MEMBERWORK/$OLCF_PROJECT/comet_work

if [ $debug = 1 ] ; then
  executable_double=$INSTALLS_DIR/install_test_$host/bin/genomics_metric
  executable_single=$INSTALLS_DIR/install_single_test_$host/bin/genomics_metric
else
  executable_double=$INSTALLS_DIR/install_release_$host/bin/genomics_metric
  executable_single=$INSTALLS_DIR/install_single_release_$host/bin/genomics_metric
fi

if [ "$metric_type" != czekanowski ] ; then
  if [ "$single" = 1 ] ; then
    executable=$executable_single
  else
    executable=$executable_double
  fi
  [[ "$sparse" == yes ]] && tag=${metric_type}_sparse || tag=${metric_type}_nonsparse
elif [ "$single" = 1 ] ; then
  executable=$executable_single
  tag=czek_single
else
  executable=$executable_double
  tag=czek_double
fi
[[ $tc != 0 ]] && tag=${tag}_tc

uid=${tag}_$(echo $(( 100000 + $num_proc_field )) | sed 's/.//')
uid=${uid}_$(echo $(( 100000 + $num_proc_vector )) | sed 's/.//')
uid=${uid}_$(echo $(( 100000 + $num_proc_repl )) | sed 's/.//')
uid=${uid}_$(echo $(( 100000 + $num_proc )) | sed 's/.//')
uid=${uid}_${num_field_local}_${num_vector}_${phase_min}_${phase_max}_${num_phase}_$$

# Output file stub
out_stub=out_3way_$uid
outfile=${out_stub}_log.txt

#ar_opts="PAMI_IBV_ENABLE_DCT=1 PAMI_ENABLE_STRIPING=1 PAMI_IBV_ADAPTER_AFFINITY=0 PAMI_IBV_QP_SERVICE_LEVEL=8 PAMI_IBV_ENABLE_OOO_AR=1"
ar_opts="PAMI_IBV_DEVICE_NAME=mlx5_0:1,mlx5_3:1 PAMI_IBV_DEVICE_NAME_1=mlx5_3:1,mlx5_0:1 PAMI_IBV_ADAPTER_AFFINITY=1 PAMI_ENABLE_STRIPING=1 PAMI_IBV_ENABLE_OOO_AR=1 PAMI_IBV_QP_SERVICE_LEVEL=8 PAMI_IBV_ENABLE_DCT=1"
launch_command="env OMP_NUM_THREADS=7 $ar_opts jsrun --smpiargs=-gpu --nrs $(( $num_node_launch * $ranks_per_node )) --bind packed:7 --cpu_per_rs 7 --gpu_per_rs 1 --rs_per_host $ranks_per_node --tasks_per_rs 1 -X 1"

#------------------------------------------------------------------------------

[[ $num_node_solve == $num_node_launch ]] &&  fastnodes_arg="" || fastnodes_arg="--fastnodes"

[[ -z "$threshold" ]] && threshold=.6
[[ -z "$metrics_shrink" ]] && metrics_shrink=50

# Command to execute, with options

if [ $metric_type != czekanowski ] ; then
  exec_command="$launch_command $executable \
    --num_way 3 \
    --num_field_local $num_field_local \
    --num_vector_local $num_vector_local \
    --metric_type $metric_type \
    --sparse $sparse \
    --all2all yes \
    --compute_method GPU \
    --problem_type $problem_type \
    --checksum $cksum \
    --num_proc_vector $num_proc_vector \
    --num_proc_field $num_proc_field \
    --num_proc_repl $num_proc_repl \
    --num_phase $num_phase --phase_min $phase_min --phase_max $phase_max \
    --num_stage $num_stage --stage_min $(( $num_stage - 1 )) \
    --threshold $threshold \
    --verbosity 1 $fastnodes_arg \
    --metrics_shrink $metrics_shrink \
    --tc $tc --num_tc_steps $num_tc_steps "
else
  exec_command="$launch_command $executable \
    --num_way 3 \
    --num_field_local $num_field_local \
    --num_vector_local $num_vector_local \
    --metric_type czekanowski \
    --all2all yes \
    --compute_method GPU \
    --problem_type $problem_type \
    --checksum $cksum \
    --num_proc_vector $num_proc_vector \
    --num_proc_field $num_proc_field \
    --num_proc_repl $num_proc_repl \
    --num_phase $num_phase --phase_min $phase_min --phase_max $phase_max \
    --num_stage $num_stage --stage_min $(( $num_stage - 1 )) \
    --verbosity 1 $fastnodes_arg"
fi

# Perform run

date | tee -a $outfile
echo "$exec_command" | tee -a $outfile
time $exec_command 2>&1 | tee -a $outfile
date | tee -a $outfile

exit

#==============================================================================
