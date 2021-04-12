#!/bin/bash

#SBATCH --account=GRONOR
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
#SBATCH --output=Network-Test.%j.out
#SBATCH --error=Network-Test.%j.err

cd /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/

# Run get info
#ech -e "\n\nGetting process information"

#echo -e "\n\nnumactl"
#srun -n 1 numactl -H

#echo -e

#echo -e "\n\ncheck-mpi cores"
#srun -n 1 which check-mpi.gnu.cori
#srun -n 8 -G 8 --cpu-bind=cores check-mpi.gnu.cori

#echo -e "\n\nhwloc-ls physical output:"
#srun -n 1 hwloc-ls -p

#echo -e "\n\nhwloc-ls logical output:"
#srun -n 1 hwloc-ls -l

#echo -e "\n\nnvidia-smi:"
#srun -n 1 nvidia-smi

export OMP_PROC_BIND=true

# Run put_to_first_core/omp_id tests
echo "Running put_to_first_core/omp_id test 1"
export OMP_NUM_THREADS=12
srun -n 1 -c 12 -G 1 --cpu-bind=verbose,mask_ldoms:0xc,0x3,0xc0,0x30 ./put_to_first_core.sh ./omp_id
echo "Done running put_to_first_core/omp_id test 1"

echo "Running put_to_first_core/omp_id test 1 24"
export OMP_NUM_THREADS=24
srun -n 1 -c 24 -G 1 --cpu-bind=verbose,mask_ldoms:0xc,0x3,0xc0,0x30 ./put_to_first_core.sh ./omp_id
echo "Done running put_to_first_core/omp_id test 1 24"

echo "Running put_to_first_core/omp_id test 4"
export OMP_NUM_THREADS=12
srun -n 4 -c 12 -G 4 --cpu-bind=verbose,mask_ldoms:0xc,0x3,0xc0,0x30 ./put_to_first_core.sh ./omp_id
echo "Done running put_to_first_core/omp_id test 4"

echo "Running put_to_first_core get_info test 4"
export OMP_NUM_THREADS=12
srun -n 4 -c 12 -G 4 --cpu-bind=verbose,mask_ldoms:0xc,0x3,0xc0,0x30 ./put_to_first_core.sh ./get_info.exe
echo "Done running put_to_first_core get_info test 4"

echo "Running put_to_first_core get_info test 4 24"
export OMP_NUM_THREADS=24
srun -n 4 -c 24 -G 4 --cpu-bind=verbose,mask_ldoms:0xc,0x3,0xc0,0x30 ./put_to_first_core.sh ./get_info.exe
echo "Done running put_to_first_core get_info test 4 24"

#export OMP_NUM_THREADS=24
#export OMP_PROC_BIND=spread
#export OMP_PLACES=cores
#export OMP_PLACES=sockets
#export OMP_PLACES=threads

# Improved settings
#echo -e "\n\n\n\nImproved settings:"
#export OMP_NUM_THREADS=24
#time srun -n 4 -c 24 -G 4 --cpu-bind=mask_ldoms:0xc,0x3,0xc0,0x30 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Improved settings 12
#echo -e "\n\n\n\nImproved settings 12:"
#export OMP_NUM_THREADS=12
#time srun -n 4 -c 12 -G 4 --cpu-bind=mask_ldoms:0xc,0x3,0xc0,0x30 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Split mask settings
#echo -e "\n\n\n\nSplit mask settings:"
#export OMP_NUM_THREADS=24
#time srun -n 4 -c 24 -G 4 --cpu-bind=mask_ldoms:0xc,0x3,0xc0,0x30 bash split_mask.sh /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Split mask settings 2
#echo -e "\n\n\n\nSplit mask settings 2:"
#export OMP_NUM_THREADS=24
#time srun -n 4 -c 24 -G 4 bash split_mask.sh /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Default settings
#echo -e "\n\nDefault settings 12:"
#export OMP_NUM_THREADS=12
#time srun -n 4 -G 4 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 12 -G 4 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Default settings
#echo -e "\n\nDefault settings 24:"
#export OMP_NUM_THREADS=24
#time srun -n 4 -c 24 -G 4 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Set cpu-binding cores
#echo -e "\n\nCPU-binding cores:"
#time srun -n 4 -G 4 --cpu-bind=cores /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 12 -G 4 --cpu-bind=cores /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 24 -G 4 --cpu-bind=cores /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Set cpu-binding sockets
#echo -e "\n\n\n\nCPU-binding sockets 12:"
#export OMP_NUM_THREADS=12
#time srun -n 4 -c 12 -G 4 --cpu-bind=sockets /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Set cpu-binding threads
#echo -e "\n\nCPU-binding sockets 24:"
#export OMP_NUM_THREADS=24
#time srun -n 4 -c 24 -G 4 --cpu-bind=sockets /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Set cpu- and gpu-binding
#echo -e "\n\nOther CPU-binding sockets 24:"
#export OMP_NUM_THREADS=24
#time srun -n 4 -c 1 -G 4 --cpu-bind=map_cpu:0-1,2-3,4-5,6-7 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

#echo -e "\n\nCPU- and GPU-binding:"
#time srun -n 4 -G 4 --cpu-bind=cores --gpu-bind=map_gpu:0,1,2,3 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 12 -G 4 --cpu-bind=cores --gpu-bind=map_gpu:0,1,2,3 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 24 -G 4 --cpu-bind=cores --gpu-bind=map_gpu:0,1,2,3 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

# Use bad gpu-binding
#echo -e "\n\nBad GPU-binding:"
#time srun -n 4 -G 4 --cpu-bind=sockets --gpu-bind=map_gpu:2,3,0,1 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 12 -G 4 --cpu-bind=sockets --gpu-bind=map_gpu:2,3,0,1 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe
#time srun -n 4 -c 24 -G 4 --cpu-bind=sockets --gpu-bind=map_gpu:2,3,0,1 /p/project/gronor/eller1/comet_dev/genomics_gpu/tools/network_info/get_info.exe

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
#time srun -n 8 -G 8 --cpu-bind=verbose,map_cpu:0,1,2,3,4,5,6,7 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 --exclusive /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

#time srun -n 128 -G 8 --cpu-bind=verbose,sockets --cpus-per-gpu 16 --gpu-bind=map_gpu:0*16,1*16,2*16,3*16,4*16,5*16,6*16,7*16 --exclusive /global/homes/p/peller/comet_dev/run_data/test_network_info/get_info.exe

