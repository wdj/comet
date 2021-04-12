#!/usr/bin/env bash
## Add mask and launch process; masks added: one GPU-aware core, the remaining cores of two NUMA domains.
## Arguments to script:
##  * None, >2: Print masks
##  * 1: Launch same process with different masks
##  * 2: Launch first process with other mask, launch second process with GPU mask
## NOTE: Highly tuned to AMD EPYC topology of JUWELS Booster. Will work for other topologies, but needs fine-tuning
## -Andreas Herten, 12 December 2020

# Aim:
# Rank 0 - Numa 2 + 3 - cores 12-23, 60-71 - GPU = 18
# Rank 1 - Numa 0 + 1 - cores  0-11, 48-59 - GPU =  6
# Rank 2 - Numa 6 + 7 - cores 36-47, 84-95 - GPU = 42
# Rank 3 - Numa 4 + 5 - cores 24-35, 72-83 - GPU = 30

## Get odd NUMA domain (which is always close to a GPU)
numa_domains=$(numactl -s | grep nodebind | sed "s/nodebind: //")
echo "numa_domains = ${numa_domains}"

for domain in $numa_domains; do
	if [ $((domain%2)) == 1 ]; then
		odd_domain=$domain
	fi
done
echo "odd_domain = ${domain}"

## GPU for odd NUMA domain
gpu_id=$(bash get_close_gpu.sh $odd_domain)
echo "gpu_id = ${gpu_id}"

## Get cores of NUMA domain
numa_cores=$(numactl -s | grep physcpubind | sed "s/physcpubind: //")
numa_cores_array=($numa_cores) ## convert to array
echo "numa_cores = ${numa_cores}"
echo "numa_cores_array = ${numa_cores_array}"

## Split list of cores into single GPU-close core one the remaining ones
N_CORES_PER_DOMAIN=6  ## this can probably be retrieved from the systems somewhere^TM
core_gpu=${numa_cores_array[$N_CORES_PER_DOMAIN]}  ## this implicitly assumes even NUMA domain before odd NUMA #domain
#core_rest=${numa_cores_array[@]/$core_gpu}
core_rest=${numa_cores_array[@]:6:6}
core_rest+=" "
core_rest+=${numa_cores_array[@]:0:6}
core_rest+=" "
core_rest+=${numa_cores_array[@]:18:6}
core_rest+=" "
core_rest+=${numa_cores_array[@]:12:6}
#core_test=("${numa_cores_array[0]}" "${numa_cores_array[1]}")
#echo "core_test = ${core_test}"
#core_rest=("${core_rest1[@]}" "${core_rest2[@]}")
core_rest_commasep=$(echo $core_rest |  tr " " ",")
echo "core_gpu = ${core_gpu}"
echo "core_rest = ${core_rest}"
echo "core_rest_commasep = ${core_rest_commasep}"

## Mask program calls to match the domains, including GPU
echo "Args = $1"
#env -u CUDA_VISIBLE_DEVICES numactl --physcpubind=$core_rest_commasep $1
CUDA_VISIBLE_DEVICES=$gpu_id numactl --physcpubind=$core_rest_commasep $1
#CUDA_VISIBLE_DEVICES=$gpu_id numactl --physcpubind=$core_gpu $1

