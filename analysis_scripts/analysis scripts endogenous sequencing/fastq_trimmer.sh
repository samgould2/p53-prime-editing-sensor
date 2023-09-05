#!/bin/bash

#cd into the correct directory
cd /Volumes/sanchezrivera/samgould/230801San/
#just run it in this folder^; doesn't seem to work otherwise...

#access the config file
#config=./config.txt
config=/Users/samgould/Documents/GitHub/p53-sensor-NGS-analysis/config_singular.txt

# Extract R1_FILE name for the current $SLURM_ARRAY_TASK_ID
i=11 #modified to start at 11 to continue running over-night (change back to 1 for full run)
max=49
while [ $i -lt $max ]
do

    R1_FILE=$(awk -v ArrayTaskID=$i '$1==ArrayTaskID {print $2}' $config)
    R2_FILE=$(awk -v ArrayTaskID=$i '$1==ArrayTaskID {print $3}' $config)
    folder_name=$(awk -v ArrayTaskID=$i '$1==ArrayTaskID {print $4}' $config)

    echo "output: $i"
    echo "R1_File: ${R1_FILE}"
    echo "R2_File: ${R2_FILE}"

    python3 /Users/samgould/Documents/GitHub/p53-sensor-NGS-analysis/fastq_trimmer.py ${R1_FILE} ${R2_FILE}

    true $(( i++ ))
done