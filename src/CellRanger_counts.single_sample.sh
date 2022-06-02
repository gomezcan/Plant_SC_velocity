#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:00:00          	# Limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1            		  # Number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1          		  # Number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=20     		# Number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=100G       	       	# Specify the real memory required per node. (20G)
#SBATCH --job-name CR_count 	    # You can give your job a name for easier identification (same as -J)


#########################################
########## 	Load Modules 	#########
module purge
#########################################

# Genome index
Index_TAIR11_CellRange=/mnt/home/gomezcan/Projects/SingleCell/Genomes/Index_TAIR11_CellRange/

# Sample name base
Sample=$1

# list of directories: multiples run/lanes from the same library
Sample_dir=$(ls -m -d ${Sample}L* | tr ' ' ',' | sed 's/,,/,/g');

echo $Sample_dir;

cellranger count --id=Counts_$Sample --transcriptome=$Index_TAIR11_CellRange \
		--fastqs=$Sample_dir \
		--sample=$Sample_dir \
        	--localcores=100 \
		--localmem=100 \
		--include-introns true
