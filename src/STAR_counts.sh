#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=0:30:00                  # Limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                       # Number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                      # Number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=100             # Number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=50G                       # Specify the real memory required per node. (20G)
#SBATCH --job-name STAR_index           # You can give your job a name for easier identification (same as -J)


########## Load Modules #########
module purge
ml GCC/10.3.0 STAR/2.7.9a
################################

for i in SRR*R1_001.fastq.gz; do
        # SRR12046049_S1_L001_R1_001.fastq.gz
        File1=$i;
        File2=${i//_R1_/_R2_};

        name=$(echo $i | tr '_' '\t' | cut -f1);
        STAR --soloType Droplet --runThreadN 100 --soloCBwhitelist 3M-february-2018.txt.gz \
                                --genomeDir ../Index_TAIR11_STAR --readFilesIn $File2 $File1 --outFileNamePrefix ${name}_star \
                                --outSAMtype BAM --soloUMIlen 12 --soloOutFileNames StarCounts_${name}
done;
