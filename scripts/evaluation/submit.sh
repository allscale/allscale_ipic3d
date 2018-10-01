#!/bin/bash

for i in 1000000 2000000 4000000 8000000 16000000 32000000 48000000 100000000
do
    #submit single node scaling experiment
    sbatch allscalecc_single_node_${i}.sh

    #submitting to normal queue for multi node weak scaling
    for n in 1 2 4 8 16 32 64
    do
        sbatch -n${n} -N${n} allscalecc_weak_${i}.sh
    done
    #submiting to big queue for multi node weak scaling
    for j in 128 256
    do
        sbatch -n$j -N$j -p big allscalecc_weak_$i.sh
    done
done
