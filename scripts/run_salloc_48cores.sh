#!/bin/bash -l

# The name of the script is myjob
#SBATCH -J allscale

#SBATCH -A 2018-24

# 100 hours wall-clock time will be given to this job
#SBATCH -t 100:00:00

# Get an email once it is ready
##SBATCH --mail-type=ALL 

#SBATCH --mem=1000000

# Number of nodes (computers)
#SBATCH --nodes=1

sleep 10d

