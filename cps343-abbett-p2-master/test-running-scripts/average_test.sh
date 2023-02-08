#!/bin/bash

# Script to run multiple trials of Power Method project code and compute
# average computation time.
#
# Assumptions
#  1. Power method project code records wall-clock time necessary to
#     estimate eigenvalue by power method.
#  2. Computation time is displayed in a line that contains the string
#     "compute".
#  3. No other output line contains the string "compute".
#  4. The computation time is the 5th white-space delimited string in
#     the line.
#  NOTE: If assumptions 2 and 4 are not correct, change the values of
#  "keyword" and "position" below.
#
# 2020-02-25 Jonathan Senning <jonathan.senning@gordon.edu>
#

if (($# < 3))
then
    echo "usage: $0 NUM_TRIALS PROGRAM MATRIX_FILE NUM_THREADS"
    exit 0
fi

# Assign parameters
num_trials=${1}
program=${2}
matrix_file=${3}
num_threads=${4:-$(nproc)}

# Change these if they are not correct for program
keyword="compute"   # this string must appear on a single line of output
position="5"        # position of time string in line

# Run tests and compute average time
AVG_TIME=$(
    for ((i=0; i<num_trials; i++))
    do
    OMP_NUM_THREADS=$num_threads $program $matrix_file |\
	grep $keyword | awk -v pos=$position '{print $pos}'
    done | awk 'BEGIN{S=0;N=0} {S=S+$1;N=N+1} END{print S/N}'
	)

# Display result
echo -n "$AVG_TIME, $num_threads,"
