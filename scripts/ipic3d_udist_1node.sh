#!/bin/bash

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2018-3-295

# The name of the script is myjob
#SBATCH -J Allscale.ipic3d.only.shared

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 00:30:00

# Number of nodes
#SBATCH --nodes 1
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=1
# Architecture
#SBATCH -C Haswell

#SBATCH -e 180924_allscale_ipic3d_udist_1node.e
#SBATCH -o 180924_allscale_ipic3d_udist_1node.o

APP=/cfs/klemming/nobackup/r/riakymch/workspace/allscale_ipic3d/build/ipic3d_allscalecc
HOME=/cfs/klemming/nobackup/r/riakymch/workspace/allscale_ipic3d

# load modules
module swap PrgEnv-cray PrgEnv-gnu
module swap gcc gcc/7.3.0
module load cmake/3.7.1 
module load git

export LD_LIBRARY_PATH=/cfs/klemming/nobackup/p/philgs/allscale/libs/ncurses-5.9/lib:$LD_LIBRARY_PATH
export PATH=/cfs/klemming/nobackup/p/philgs/allscale/libs/ruby-1.9.3-p125/bin:$PATH
export CRAYPE_LINK_TYPE=dynamic
export CRAY_ROOTFS=DSL

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cfs/klemming/nobackup/r/riakymch/allscale_compiler/build/allscale_runtime-prefix/src/allscale_runtime-build/src:/cfs/klemming/nobackup/r/riakymch/allscale_compiler/build/third_party/boost/lib:/cfs/klemming/nobackup/r/riakymch/allscale_compiler/build/hpx-prefix/src/hpx-build/lib


export ALLSCALE_MONITOR=0
export ALLSCALE_RESILIENCE=0
echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"
PARTICLES=$((8000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((16000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((32000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((48000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32


export ALLSCALE_MONITOR=1
export ALLSCALE_RESILIENCE=0
echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"
PARTICLES=$((8000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((16000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((32000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((48000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32


export ALLSCALE_MONITOR=0
export ALLSCALE_RESILIENCE=1
echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"
PARTICLES=$((8000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((16000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((32000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((48000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32


export ALLSCALE_MONITOR=1
export ALLSCALE_RESILIENCE=1
echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"
PARTICLES=$((8000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((16000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((32000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
PARTICLES=$((48000000 * SLURM_NNODES))
aprun -n 1 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32
