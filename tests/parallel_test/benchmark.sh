#!/bin/bash

# Name of your executable
EXE="../../roda.exe"

INPUT_FILE="input.inp"
# File to save results
OUTPUT="benchmark_results.txt"

echo "Threads | Time (s)" > $OUTPUT
echo "-------------------" >> $OUTPUT

# Loop through thread counts: 1, 2, 4, 8, 12, 16
for t in 1 2 4 8 10 14 16
do
    echo "Running with $t threads..."
    
    # Set the number of threads
    export OMP_NUM_THREADS=$t
    
    # Run and capture the "Total CPU time" line from your output
    # (Assuming you updated it to print Wall Time, or we use the linux 'time' command)
    
    # We use the Linux 'time' command as a backup to measure real time
    # %e gives real time in seconds
    /usr/bin/time -f "%e" $EXE $INPUT_FILE > temp_output.log 2> time_log.txt
    
    REAL_TIME=$(cat time_log.txt)
    
    echo "$t       | $REAL_TIME" >> $OUTPUT
    echo "Finished $t threads in $REAL_TIME seconds."
done

echo "Benchmark complete. Results in $OUTPUT"
cat $OUTPUT

