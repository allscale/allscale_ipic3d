#!/bin/bash

ITER1=100
ITER2=99

for N in 1000000 2000000 4000000 8000000 16000000 32000000 48000000 100000000
do
    echo $N
    for f in allscalecc_ipic3d_weak_${N}.o*
    do
        TIME_STR=$(grep "Simulation took" $f | head -n1)
        TIME_MON_STR=$(grep "Simulation took" $f | tail -n1)
        TMP=$(grep "Throughput:" $f | tail -n1 | awk -F':' '{print $2}' | awk -F' ' '{print $1}')
        TMP=$(python -c "print($TMP)")

        PARTICLES_STR=$(grep "particles / cell:" $f | head -n1 | awk -F':' '{print $2}' | sed 's/(\(.*\))/\1/')


        NNODES=$(grep "Requested resources" $f | awk -F':' '{print $2}' | awk -F',' '{print $2}' | awk -F'=' '{print $2}')
        PARTICLES="$N * $NNODES"
        #$(echo "($PARTICLES_STR) * (16 * 16 * 16) * 64" | bc)

        TIME1=$(echo $TIME_STR | awk -F' ' '{print $3}')
        TIME2=$(echo $TIME_STR | awk -F' ' '{print $4}' | sed 's/^(\(.*\))$/\1/')
        TIME3=$(echo $TIME_MON_STR | awk -F' ' '{print $3}')
        TIME4=$(echo $TIME_MON_STR | awk -F' ' '{print $4}' | sed 's/^(\(.*\))$/\1/')

        TIME1=$(echo "min($TIME1, $TIME3)" | bc min.bc)
        TIME2=$(echo "min($TIME2, $TIME4)" | bc min.bc)


        THROUGHPUT1=$(echo "($PARTICLES * $ITER1) / $TIME1" | bc -l)
        THROUGHPUT2=$(echo "($PARTICLES * $ITER2) / $TIME2" | bc -l)
        PERCENTAGE=$(echo "scale=2; 100 * ($TIME1 - $TIME2)/$TIME1" | bc )
        PERCENTAGE2=$(echo "scale=2; 100 * ($THROUGHPUT2 - $THROUGHPUT1)/$THROUGHPUT2" | bc )
        TEST=$(echo "$THROUGHPUT1 - $TMP" | bc -l)

        echo "$NNODES, $THROUGHPUT1, $THROUGHPUT2, $PERCENTAGE, $PERCENTAGE2"
        echo "-----"
    done
done
