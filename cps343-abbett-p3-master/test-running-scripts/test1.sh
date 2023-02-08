#!/bin/bash

max_task=180
max_node=16
max_proc_per_node=12

# For the AVG_TIME code
keyword="compute"
position="5"
num_trials=5

program="./proj3"
matrix_file="/gc/cps343/matrix/A-10000x10000.dat"

echo "Test 1: num_tasks, avgtime"

singlecommand="salloc -Q --ntasks=1 --exclusive mpiexec $program $matrix_file"
    AVG_TIME=$(
               for ((i=0; i<num_trials; i++))
               do
                   $singlecommand |\
                   grep $keyword | awk -v pos=$position '{print $pos}'
               done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                )
            echo -n "1, $AVG_TIME,"
            echo

for ((num_tasks_increment=1; num_tasks_increment<=9; num_tasks_increment++ ))
do 
    num_tasks=$(($num_tasks_increment*5))
    if ((num_tasks != 45))
    then
        command="salloc -Q --ntasks=$num_tasks --exclusive mpiexec $program $matrix_file"
        AVG_TIME=$(
                for ((i=0; i<num_trials; i++))
                do
                    $command |\
                    grep $keyword | awk -v pos=$position '{print $pos}'
                done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                    )
                echo -n "$num_tasks, $AVG_TIME,"
                echo
    else
            command="salloc -Q --ntasks=44 --exclusive mpiexec $program $matrix_file"
            AVG_TIME=$(
                for ((i=0; i<num_trials; i++))
                do
                    $command |\
                    grep $keyword | awk -v pos=$position '{print $pos}'
                done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                    )
                echo -n "$num_tasks, $AVG_TIME,"
                echo
    fi
done