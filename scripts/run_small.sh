#!/bin/bash -l
# The -l above is required to get the full environment with modules

export INPUT=../inputs/small.inp
export CMD='/usr/bin/time --format="%e" ./../build/app/ipic3d'

# Run the executable named myexe 
export NUM_WORKERS=48
$CMD $INPUT > 1.txt

export NUM_WORKERS=32
$CMD $INPUT > 1.txt

export NUM_WORKERS=16
$CMD $INPUT > 1.txt

export NUM_WORKERS=8
$CMD $INPUT > 1.txt

export NUM_WORKERS=4
$CMD $INPUT > 1.txt

export NUM_WORKERS=2
$CMD $INPUT > 1.txt

export NUM_WORKERS=1
$CMD $INPUT > 1.txt
