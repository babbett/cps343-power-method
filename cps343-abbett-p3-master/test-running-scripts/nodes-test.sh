#!/bin/bash

# For the AVG_TIME code
keyword="compute"
position="5"
num_trials=5

program=$1
matrix_file="/gc/cps343/matrix/A-10000x10000.dat"

echo "$program num_nodes, tasks_per, avgtime"
for ((num_nodes=1; num_nodes<=11; num_nodes++ ))
do 
    for ((tasks_per=1; tasks_per<=4; tasks_per++))
    do
	    command="salloc -Q --nodes=$num_nodes --ntasks-per-node=$tasks_per --exclusive mpiexec $program $matrix_file"
	    AVG_TIME=$(
		       for ((i=0; i<num_trials; i++))
		       do
			   $command |\
			   grep $keyword | awk -v pos=$position '{print $pos}'
		       done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
			)
		    echo "$num_nodes $tasks_per $AVG_TIME"
    done
done