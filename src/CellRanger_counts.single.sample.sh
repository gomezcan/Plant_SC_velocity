#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:00:00          	# Limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1            		  # Number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1          		  # Number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=20     		# Number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=100G       	       	# Specify the real memory required per node. (20G)
#SBATCH --job-name CR_count 	    # You can give your job a name for easier identification (same as -J)


########## Load Modules #########


module purge

File=/mnt/home/gomezcan/Projects/SingleCell/RawData

cellranger count --id=$1 \
        --transcriptome=/mnt/home/gomezcan/Projects/SingleCell/Genomes/Index_TAIR11_CellRange/ \
        --fastqs=$File/$1/ \
        --sample=$1 \
	--include-introns
        --localcores=20 --localmem=100
