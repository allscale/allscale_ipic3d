#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2018-3-295

# The name of the script is myjob
#SBATCH -J Allscale.ipic3d.shared

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 01:00:00

# Number of nodes
#SBATCH -N 1
# Number of MPI processes.
#SBATCH -n 1
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=32
# Architecture
#SBATCH -C Haswell

#SBATCH -e error_file.e
#SBATCH -o output_file.o

# Run the executable named myexe 
# and write the output into my_output_file
