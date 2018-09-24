#!/bin/bash

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2018-3-295

# The name of the script is myjob
#SBATCH -J Allscale.ipic3d.only.shared

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 04:00:00

# Number of nodes
#SBATCH --nodes 128
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=1
# Architecture
#SBATCH -C Haswell

#SBATCH -e 180924_allscale_ipic3d_udist_128nodes.e
#SBATCH -o 180924_allscale_ipic3d_udist_128nodes.o

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


for RES in 0 1
do
	for MON in 0 1
	do
		export ALLSCALE_MONITOR=$MON
		export ALLSCALE_RESILIENCE=$RES

        	echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"

		for PARTICLES in 8000000 16000000 32000000 48000000
		do
			PARTICLES=$((PARTICLES * SLURM_NNODES))
			aprun -n 128 -N 1 -d 32 $APP :U:$PARTICLES --hpx:threads=32 
		done

	done
done
