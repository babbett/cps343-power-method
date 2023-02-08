#!/bin/bash

max_task=180
max_node=16
max_proc_per_node=12

# For the AVG_TIME code
keyword="compute"
position="5"
num_trials=5

program="./proj3_omp"
matrix_file="/gc/cps343/matrix/A-10000x10000.dat"

echo "Test 1: num_tasks, cpus_per, avgtime"
for ((num_tasks=1; num_tasks<=44; num_tasks++ ))
do 
    for ((num_cpus=1; num_cpus<=4; num_cpus++))
    do
	    if (( $(($num_cpus*$num_tasks)) <= 44))
	    then
		    command="salloc -Q --ntasks=$num_tasks --cpus-per-task=$num_cpus --exclusive mpiexec $program $matrix_file"
		    AVG_TIME=$(
			       for ((i=0; i<num_trials; i++))
			       do
				   $command |\
				   grep $keyword | awk -v pos=$position '{print $pos}'
			       done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
				)
			    echo -n "$num_tasks, $num_cpus, $AVG_TIME,"
			    echo
	    fi
    done
done