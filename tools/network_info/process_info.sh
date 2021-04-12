#!/usr/bin/env bash
echo -n "MPI Rank: $MPI_LOCALRANKID;"
echo -n "CUDA_VIS_DEV: $CUDA_VISIBLE_DEVICES;"
echo -n "$(taskset -c -p $$);"
echo -n "Last CPU core: $(awk '{print $39}' /proc/$$/stat)"

echo ""

