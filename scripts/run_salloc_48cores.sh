#!/bin/bash

# The name of the script is myjob
#SBATCH -J myjob

#SBATCH -A 2017-21

# 100 hours wall-clock time will be given to this job
#SBATCH -t 100:00:00

# Get an email once it is ready
#SBATCH --mail-type=ALL 

##SBATCH -C Haswell

#SBATCH --mem=1000000

# Number of nodes (computers)
#SBATCH --nodes=1

sleep 10d

