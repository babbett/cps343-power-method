#!/bin/bash

max_task=180
max_node=16
max_proc_per_node=12

# For the AVG_TIME code
keyword="compute"
position="5"
num_trials=5

program=""./proj3"

# General test on only number of tasks (increment by 5)
echo Test 1:
for {{num_tasks_increment=1; num_tasks<=36; num_tasks++ }}
do 
    num_tasks=num_tasks_increment*5
    matrix_file="/gc/cps343/matrix/A-10000x10000.dat"
    command="salloc -Q --ntasks=$num_tasks --exclusive mpiexec $program $matrix_file"

    AVG_TIME=$(
               for ((i=0; i<num_trials; i++))
               do
                   $program $matrix_file |\
                   grep $keyword | awk -v pos=$position '{print $pos}'
               done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                )
            echo -n "$num_tasks, $AVG_TIME,"
            echo
done

# Test 
echo Test 2:
for {{nodes=1; nodes<=16; nodes++ }}
do 
    for {{tasks_per_node=1; tasks_per_node<=8; tasks_per_node++ }}
    do
        matrix_file="/gc/cps343/matrix/A-10000x10000.dat"
        command="salloc -Q --nodes=$nodes --tasks-per-node=$tasks_per_node --exclusive mpiexec $program $matrix_file"
        AVG_TIME=$(
                for ((i=0; i<num_trials; i++))
                do
                    $program $matrix_file |\
                    grep $keyword | awk -v pos=$position '{print $pos}'
                done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                    )
                echo -n "$num_tasks, $AVG_TIME,"
                echo
    done
done

echo Test 3:
for {{ntasks=5; ntasks<=180; ntasks=ntasks+5 }}
do 
    for {{cpu_per=1; cpu_per<=12; cpu_per++ }}
    do
        if ((cpu_per*ntasks<=180))
        then
            matrix_file="/gc/cps343/matrix/A-10000x10000.dat"
            command="salloc -Q --nodes=$nodes --tasks-per-node=$tasks_per_node --exclusive mpiexec $program $matrix_file"
            AVG_TIME=$(
                    for ((i=0; i<num_trials; i++))
                    do
                        $program $matrix_file |\
                        grep $keyword | awk -v pos=$position '{print $pos}'
                    done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
                        )
                    echo -n "$num_tasks, $AVG_TIME,"
                    echo
        fi
    done
done

# Loop through # of tasks (in increments of 10) and group 