#!/bin/bash

# Script to run average_test trials over different matrix sizes and programs
#
# Created by Ben Abbett
num_trials=5
keyword="compute"
position="5"
maxproc=16

for program in "./proj3" "./proj3_omp"
do
    echo "On program: $program "
    for size in 50 100 500 1000 2000 5000 10000
    do
        
        for ((nodes=1; nodes < maxproc; nodes++))
        do 
            matrix_file="/gc/cps343/matrix/A-${size}x${size}.dat"
            command="salloc -Q --nodes=$nodes mpiexec $program $matrix_file"
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

# for program in "./proj1_forloop" "./proj1_cblas" "./proj2_forloop_simple" 
# do
#     echo "On program: $program "
#     for size in 5 10 50 100 500 1000 2000 5000 10000 
#     do
#         matrix_file="/gc/cps343/matrix/A-${size}x${size}.dat"
#         AVG_TIME=$(
#                 for ((i=0; i<num_trials; i++))
#                 do
#                 $program $matrix_file |\
#                 grep $keyword | awk -v pos=$position '{print $pos}'
#                 done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
#                 )
#         echo -n "$AVG_TIME, $size, ,"
#         echo
#     done
# done

# for program in "./proj2_forloop" "./proj2_cblas"
# do
#     echo "On program: $program "
#     for size in 5 10 50 100 500 1000 2000 5000 10000 
#     do
#         matrix_file="/gc/cps343/matrix/A-${size}x${size}.dat"
#         for ((num_threads=1; num_threads<5; num_threads++))
#         do
#             ./average_test.sh $num_trials $program $matrix_file $num_threads
#             echo " $size,"
#         done
#     done
# done