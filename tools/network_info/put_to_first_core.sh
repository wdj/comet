#!/usr/bin/env bash

## Get odd NUMA domain (which is always close to a GPU)
numa_domains=$(numactl -s | grep nodebind | sed "s/nodebind: //")

for domain in $numa_domains; do
        if [ $((domain%2)) == 1 ]; then
                odd_domain=$domain
        fi
done

## GPU for odd NUMA domain
gpu_id=$(bash get_close_gpu.sh $odd_domain)

## Get cores of NUMA domain
numa_cores=$(numactl -s | grep physcpubind | sed "s/physcpubind: //")
numa_cores_array=($numa_cores) ## convert to array

## Split list of cores into single GPU-close core one the remaining ones
N_CORES_PER_DOMAIN=6  ## this can probably be retrieved from the systems somewhere^TM
core_gpu=${numa_cores_array[$N_CORES_PER_DOMAIN]}  ## this implicitly assumes even NUMA domain before odd NUMA domain
core_rest=${numa_cores_array[@]/$core_gpu}

core_resorted=( $core_gpu $core_rest )

_OMP_PLACES=""
for i in "${core_resorted[@]}"
do
	_OMP_PLACES+="{$i},"
done

if [[ -z "$OMP_PLACES" ]]; then
	export OMP_PLACES=${_OMP_PLACES::-1}  # removing trailing comma
else
	echo "OMP_PLACES already set. Not touching it."
fi

#echo $OMP_PLACES

$1

