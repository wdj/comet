#!/usr/bin/env bash

declare -A gpus=( ["3"]="0" ["1"]="1" ["7"]="2" ["5"]=3)

echo "${gpus[$1]}"

