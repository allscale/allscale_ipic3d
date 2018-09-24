#!/bin/bash -l

#SBTACH -q work
#SBATCH -t 02:00:00
#SBATCH -J allscalecc_ipic3d_weak_8000000

. ~/.bashrc
module list

APP=$HOME/build/meggie/release/allscale/examples/ipic3d_allscalecc_data

BASE_PARTICLES=8000000
PARTICLES=$((BASE_PARTICLES * SLURM_NNODES))

echo "Running on $i cores, with $BASE_PARTICLES * $SLURM_NNODES = $PARTICLES particles"
srun $APP --hpx:threads=20 -Ihpx.max_background_threads=10 :U:$PARTICLES -Ihpx.stacks.use_guard_pages=0

echo "Running on $i cores, with $BASE_PARTICLES * $SLURM_NNODES = $PARTICLES particles with monitoring"
export ALLSCALE_MONITOR=1
srun $APP --hpx:threads=20 :U:$PARTICLES -Ihpx.stacks.use_guard_pages=0
