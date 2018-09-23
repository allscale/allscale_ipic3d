#!/bin/bash

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2018-3-295

# The name of the script is myjob
#SBATCH -J Allscale.ipic3d.only.shared

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 04:00:00

# Number of nodes
#SBATCH -N 1
# Number of MPI processes.
#SBATCH -n 1
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=32
# Architecture
#SBATCH -C Haswell

#SBATCH -e allscale_ipic3d_single_node01.e
#SBATCH -o allscale_ipic3d_single_node01.o

APP=/cfs/klemming/nobackup/r/riakymch/workspace/allscale_ipic3d/build/ipic3d_allscalecc
HOME=/cfs/klemming/nobackup/r/riakymch/workspace/allscale_ipic3d

#submit single node scaling experiment
export ALLSCALE_MONITOR=0
export ALLSCALE_RESILIENCE=1

#for OPTION in U E C
#do
#    for i in 1 2 4 8 16 32
#    do
#        for PARTICLES in 1000000 2000000 4000000 8000000 16000000 32000000 64000000
#        do
#            echo "Running on $i cores, with $PARTICLES particles for $OPTION distribution"
#            aprun -n 1 -N 1 $APP --hpx:threads=$i :$OPTION:$PARTICLES -Ihpx.stacks.use_guard_pages=0
#        done
#    done
#done

#for RES in 0 1
#do
#	for MON in 0 1
#	do
#		export ALLSCALE_MONITOR=$MON
#		export ALLSCALE_RESILIENCE=$RES

        echo "Running on ipic3d with MONITORING=$ALLSCALE_MONITOR and RESILIENCE=$ALLSCALE_RESILIENCE"

		# 1 core
		aprun -n 1 -N 1 $APP --hpx:threads=1 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_1c.inp
		# 2 cores
		aprun -n 1 -N 1 $APP --hpx:threads=2 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_2c.inp
		# 4 cores
		aprun -n 1 -N 1 $APP --hpx:threads=4 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_4c.inp
		# 8 cores
		aprun -n 1 -N 1 $APP --hpx:threads=8 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_8c.inp
		# 16 cores
		aprun -n 1 -N 1 $APP --hpx:threads=16 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_16c.inp
		# 32 cores
		aprun -n 1 -N 1 $APP --hpx:threads=32 -Ihpx.stacks.use_guard_pages=0 $HOME/inputs/shared_32c.inp

#	done
#done


