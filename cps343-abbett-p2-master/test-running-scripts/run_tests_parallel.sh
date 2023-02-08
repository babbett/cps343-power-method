#!/bin/bash

# Run a tests on a specific parallel program for n iterations

num_trials=${1}
program=${2}

for size in 5 10 50 100 500 1000 2000 5000 10000 
do
    matrix_file="/gc/cps343/matrix/A-${size}x${size}.dat"
    for ((num_threads=1; num_threads<5; num_threads++))
    do
        ./average_test.sh $num_trials $program $matrix_file $num_threads
        echo " $size,"
    done
done