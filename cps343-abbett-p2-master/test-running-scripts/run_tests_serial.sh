#!/bin/bash

# Run a tests on a specific program for n iterations

num_trials=${1}
program=${2}

keyword="compute"   # this string must appear on a single line of output
position="5"        # position of time string in line       

for size in 5 10 50 100 500 1000 2000 5000 10000 
do
    matrix_file="/gc/cps343/matrix/A-${size}x${size}.dat"
    AVG_TIME=$(
            for ((i=0; i<num_trials; i++))
            do
            $program $matrix_file |\
            grep $keyword | awk -v pos=$position '{print $pos}'
            done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
            )
    echo -n "$AVG_TIME, $size, ,"
    echo
done