#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2017-21

# The name of the script is myjob
#SBATCH -J AllScale

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 04:00:00

# Number of nodes
#SBATCH -N 1

# Nodes with 4 NUMA domains
#SBATCH --mem=1000000

#SBATCH -e error_tine_0607.o
#SBATCH -o output_tiny_0607.o

module load gcc/5.3.0

#export INPUT=../inputs/testMagnetosphere3Dsmall_orig.inp
export INPUT=../inputs/tiny.inp
export CMD="/usr/bin/time ./../build/app/ipic3d"
#export CMD="/usr/bin/time --format='%C %e' ./../build/app/ipic3d"

# Run the executable named myexe 
export NUM_WORKERS=32
$CMD $INPUT

export NUM_WORKERS=16
$CMD $INPUT

export NUM_WORKERS=8
$CMD $INPUT

export NUM_WORKERS=4
$CMD $INPUT

export NUM_WORKERS=2
$CMD $INPUT

export NUM_WORKERS=1
$CMD $INPUT
