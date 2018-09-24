#!/bin/bash -l

#SBATCH -n1 -N1
#SBTACH -q work
#SBATCH -t 02:00:00
#SBATCH -J allscalecc_ipic3d_single_node

. ~/.bashrc
module list

APP=$HOME/build/meggie/release/allscale/examples/ipic3d_allscalecc_data

PARTICLES=1000000

for i in 1 2 4 8 10 20
do
    echo "Running on $i cores, with $PARTICLES particles"
    srun $APP --hpx:threads=$i :U:$PARTICLES -Ihpx.stacks.use_guard_pages=0
done
