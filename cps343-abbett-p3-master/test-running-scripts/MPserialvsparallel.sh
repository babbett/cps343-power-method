#!/bin/bash

# Script to run average_test trials over different matrix sizes and programs
#
# Created by Ben Abbett
num_trials=5
keyword="compute"
position="5"
maxproc=10

for program in "./proj3" "./proj3_omp"
do
    echo "On program: $program "
    for size in 50 100 500 1000 2000 5000 10000
    do
        
        for ((nodes=1; nodes < maxproc; nodes++))
        do 
            matrix_file="/gc/cps343/matrix/A-${size}x${size}.dat"
            command="salloc -Q --exclusive --nodes=$nodes mpiexec $program $matrix_file"
            AVG_TIME=$(
                for ((i=0; i<num_trials; i++))
                do
                    $program $matrix_file |\
                    grep $keyword | awk -v pos=$position '{print $pos}'
                done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                )
            echo -n "$nodes, $size, $AVG_TIME"
            echo
        done
    done
done