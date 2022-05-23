#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=5:00:00                  # Limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                       # Number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=100                    # Number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1               # Number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=20G                       # Specify the real memory required per node. (20G)
#SBATCH --job-name sra00                # You can give your job a name for easier identification (same as -J)

########## Load Modules #########

module purge
ml GCCcore/8.2.0 parallel-fastq-dump/0.6.5-Python-3.7.2

### $1 tab limited file, Example:
#SRR12046122    Paired
#SRR12046123    Paired
#SRR12046124    Paired

while read -r -a line; do
        if [[ ${line[1]} == "Single" ]];  then

                parallel-fastq-dump -s ${line[0]} -t 100 --gzip -O ./ ;

        elif [[ ${line[1]} == "Paired" ]]; then
                parallel-fastq-dump -s ${line[0]} -t 100 --split-files --gzip -O ./;


        fi

done < $1
