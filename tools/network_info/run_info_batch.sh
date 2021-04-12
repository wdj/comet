#!/bin/bash

#SBATCH --constraint=dgx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 128
#SBATCH --gpus=8
#SBATCH --cpus-per-task=1
#SBATCH --time=15
#SBATCH --account=m1759

# Run get info
echo -e "\n\nGetting process information"

echo -e "\n\nnumactl"
srun -n 1 numactl -H

#echo -e

#echo -e "\n\ncheck-mpi cores"
#srun -n 1 which check-mpi.gnu.cori
#srun -n 8 -G 8 --cpu-bind=cores check-mpi.gnu.cori

echo -e "\n\nhwloc-ls physical output:"
srun -n 1 hwloc-ls -p

echo -e "\n\nhwloc-ls logical output:"
srun -n 1 hwloc-ls -l

echo -e "\n\nnvidia-smi:"
srun -n 1 nvidia-smi

# Each process has access to 1 GPU, but only 4 seem to be used
#time srun -n 128 -G 8 --cpu-bind=verbose,cores --gpu-bind=closest /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Each process has access to all 8 GPUs
#time srun -n 128 -G 8 --cpu-bind=verbose,cores /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Each process is bound to a different GPU
#time srun -n 8 -G 8 --cpu-bind=verbose,cores --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Each process is bound to 1 GPU, but not all 8 are used, and the ones that are used are used in different amounts
#time srun -n 8 -G 8 --cpu-bind=verbose,cores --gpu-bind=closest /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Binds 16 ranks to each GPU (0-15,16-31,etc), but cores are different from ranks and don't line up well
#time srun -n 128 -G 8 --cpu-bind=verbose,cores --gpu-bind=map_gpu:0*16,1*16,2*16,3*16,4*16,5*16,6*16,7*16 /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Similar to previous, but cores are still different from ranks
#time srun -n 128 -G 8 --cpu-bind=verbose,rank_ldom --gpu-bind=map_gpu:0*16,1*16,2*16,3*16,4*16,5*16,6*16,7*16 /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Similar to previous, but cores are still different from ranks
#time srun -n 128 -G 8 --cpu-bind=verbose,rank --gpu-bind=map_gpu:0*16,1*16,2*16,3*16,4*16,5*16,6*16,7*16 /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

# Bound each gpu to different rank, with 4 ranks on 1 cpu and 4 ranks on 2nd
time srun -n 8 -G 8 --cpu-bind=verbose,map_cpu:0,1,2,3,4,5,6,7 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 --exclusive /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

#time srun -n 128 -G 8 --cpu-bind=verbose,sockets --cpus-per-gpu 16 --gpu-bind=map_gpu:0*16,1*16,2*16,3*16,4*16,5*16,6*16,7*16 --exclusive /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

