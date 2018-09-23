#!/bin/bash

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2018-3-295

# The name of the script is myjob
#SBATCH -J Allscale.ipic3d.only.shared

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 08:00:00

# Number of nodes
#SBATCH --nodes 4
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=32
# Architecture
#SBATCH -C Haswell

#SBATCH -e allscale_ipic3d_512K_4nodes.e
#SBATCH -o allscale_ipic3d_512K_4nodes.o

APP=/cfs/klemming/nobackup/r/riakymch/workspace/allscale_ipic3d/build/ipic3d_allscalecc
HOME=/cfs/klemming/nobackup/r/riakymch/workspace/allscale_ipic3d

#submit single node scaling experiment


for RES in 0 1
do
	for MON in 0 1
	do
		export ALLSCALE_MONITOR=$MON
		export ALLSCALE_RESILIENCE=$RES

        echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"

		aprun -n 4 -N 1 $APP --hpx:threads=32 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_8c.inp

	done
done
