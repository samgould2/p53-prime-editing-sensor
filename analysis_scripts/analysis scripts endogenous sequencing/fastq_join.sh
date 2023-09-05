#!/bin/bash

#cd into the correct directory
cd /Volumes/sanchezrivera/samgould/230801San/
#just run it in this folder^; doesn't seem to work otherwise...

#access the config file
#config=./config.txt
config=/Users/samgould/Documents/GitHub/p53-sensor-NGS-analysis/config_singular.txt

# Extract R1_FILE name for the current $SLURM_ARRAY_TASK_ID
i=1
max=49
while [ $i -lt $max ]
do

    R1_FILE=$(awk -v ArrayTaskID=$i '$1==ArrayTaskID {print $2}' $config)
    R2_FILE=$(awk -v ArrayTaskID=$i '$1==ArrayTaskID {print $3}' $config)
    folder_name=$(awk -v ArrayTaskID=$i '$1==ArrayTaskID {print $4}' $config)

    echo "output: $i"
    echo "folder: ${folder_name}"
    echo "R1_File: ${R1_FILE}"
    echo "R2_File: ${R2_FILE}"

    fastq-join ${R1_FILE} ${R2_FILE} -o ./fastq-join/${folder_name}_%.fastq

    true $(( i++ ))
done

#R1_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
#R2_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
#folder_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)


#fastq-join ${R1_FILE} ${R2_FILE} -o ./fastq_join/${folder_name}_%.fastq

#test 
#echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${R1_FILE} and ${R2_FILE} and the output folder is ${folder_name}" >> output.txt