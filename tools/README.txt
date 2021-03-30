
Genomics calculation workflow
=============================

NOTE: see https://code.ornl.gov/wjd/genomics_gpu/tree/master/tools
for the reference copy of the files in this directory

Preprocessing example
---------------------

ssh andes
export OLCF_PROJECT=stf006
salloc -N1 -A $OLCF_PROJECT -t 1440

TOOLS_DIR=$HOME/genomics/genomics_gpu/tools
pushd $TOOLS_DIR
./make.sh
popd

export PATH="${PATH}:$TOOLS_DIR"

shuf Area-0-Sector-10074.txt > Area-0-Sector-10074_shuf.txt

preprocess duo Area-0-Sector-10074_shuf.txt Area-0-Sector-10074_shuf.bin

line_labels.sh Area-0-Sector-10074_shuf.txt Area-0-Sector-10074_shuf_line_labels.txt

line_indices Area-0-Sector-10074_shuf.txt Area-0-Sector-10074_shuf_line_indices.bin

sed -e 's/.*/AT/' <Area-0-Sector-10074_shuf_line_labels.txt >Area-0-Sector-10074_shuf_allele_labels.txt

OPTIONAL:

preprocess_validate duo Area-0-Sector-10074_shuf.bin Area-0-Sector-10074_shuf.txt Area-0-Sector-10074_shuf_allele_labels.txt

How to select CoMet settings, 2-way case
----------------------------------------

When using a GPU-capable system like Titan or Summit, the CoMet code is
set up to use 1 GPU per MPI rank (denoted below as "proc" or "process").
On Titan this means 1 rank per node, on Summit 6 ranks per node.
OpenMP is also used to divide up the CPU cores on the node to
speed up the non-GPU work (16 OpenMP threads per MPI rank on Titan,
7 threads per rank on Summit).

CoMet has several tuning parameters for running in parallel:

  --num_proc_vector: the number of process along the "vector" axis.

  --num_proc_field: the number of process along the "field" axis.
    Note num_proc_vector X num_proc_field = num_proc equals the total
    number of MPI processes used to solve the problem.  num_proc
    nodes should be selected by qsub/aprun on Titan.  On Summit which
    has 6 GPUs per node, one should select ceil(num_proc/6) nodes.

  --num_phase, --phase_min, --phase_max: Since solving an entire
    problem meay require more memory than is available, it is
    possible to divide the computation into smaller "phases,"
    each of which computes only a subset of the result metrics.
    Computing a single phase fully utilizes all GPUs in the
    allocation.  Example: one could set num_phase=200,
    phase_min=0, phase_max=0 to compute the first phase only
    out of 200 phases.  One can num_phase=200, phase_min=0,
    phase_max=199 to compute all phases, while only requiring
    the amount of memory required to compute one phase at a time.
    NOTE: For computing all the metrics for a given fixed problem,
    one must fix the num_proc_vector and num_phase values for all phases
    computed to ensure a consistent definition of which metrics
    are contained in each phase.

The rationale for adjusting the settings is:
  1) to make the number of vectors and vector elements on each GPU as large
     as possible to mximize performance;
  2) to make num_phase large enough to ensure metrics will fit onto memory;
  3) to run at least several phases per individual run to amortize setup
     costs.

The spreadsheet CoMet_settings_Tool.xlsx can be used to assist with determining
good settings.


CoMet execution example, 2-way case
-----------------------------------

Summit:
ssh summit.olcf.ornl.gov
bsub -P $OLCF_PROJECT -Is -nnodes 1 -W 5 i_is $SHELL

executable=$HOME/genomics/install_single_release_summit/bin/genomics_metric ar_opts='PAMI_IBV_ENABLE_DCT=1 PAMI_ENABLE_STRIPING=1 PAMI_IBV_ADAPTER_AFFINITY=0 PAMI_IBV_QP_SERVICE_LEVEL=8 PAMI_IBV_ENABLE_OOO_AR=1'

launch_command="env OMP_NUM_THREADS=7 $ar_opts jsrun --nrs 1 --bind packed:7 --cpu_per_rs 7 --gpu_per_rs 1 --rs_per_host 1 --tasks_per_rs 1 -X 1"

$launch_command $executable --num_way 2 --metric_type duo --sparse yes --num_vector 3793 --num_field 414640 --all2all yes --compute_method GPU --num_proc_vetor 1 --tc 1 --num_tc_steps 4 --verbosity 1 --checksum no --threshold .7 --input_file Area-0-Sector-10074_shuf.bin --output_file_stub Area-0-Sector-10074_shuf-comet


Postprocessing example
----------------------

num_way=2

postprocess $num_way Area-0-Sector-10074_shuf_allele_labels.txt \
  Area-0-Sector-10074_shuf_line_labels.txt \
  Area-0-Sector-10074_shuf-comet_0.bin \
  Area-0-Sector-10074_shuf-comet_0.txt

tr ' ' '\t' < Area-0-Sector-10074-comet_shuf_0.txt | \
cut -f1-$(( 2 * $num_way )) | \
validate duo $num_way Area-0-Sector-10074_shuf.txt Area-0-Sector-10074_shuf_line_indices.bin > Area-0-Sector-10074_shuf_validate.txt

