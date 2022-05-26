#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=0:30:00                  # Limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                       # Number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                      # Number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=100             # Number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=50G                       # Specify the real memory required per node. (20G)
#SBATCH --job-name STAR_mapping         # You can give your job a name for easier identification (same as -J)


########## Load Modules #########
module purge
ml GCC/10.3.0 STAR/2.7.9a
################################

START_v2() {
        ## process 10x chemistry v3 of 3' 10x
        
        ## Files input: file R1, eg: SRR12046049_S1_L001_R1_001.fastq.gz
        File1=$1; 
        File2=${i//_R1_/_R2_};

        name=$(echo $i | tr '_' '\t' | cut -f1); # set sample name
        
        # mapping
        STAR --soloType Droplet --runThreadN 60 \
                                --soloCBwhitelist 737K-august-2016.txt \
                                --genomeDir Index_TAIR11_STAR \
                                --readFilesIn $File2 $File1 --outFileNamePrefix ${name}_star. \
                                --outSAMtype BAM Unsorted --soloOutFileNames Counts_${name} --readFilesCommand gunzip -c \
                                --soloCBstart 1 \
                                --soloCBlen 16 \
                                --soloUMIstart 17 \
                                --soloUMIlen 10 

}

START_v3() {
        # SRR12046049_S1_L001_R1_001.fastq.gz
        File1=$i;
        File2=${i//_R1_/_R2_};

        name=$(echo $i | tr '_' '\t' | cut -f1);
        STAR --soloType Droplet --runThreadN 100 --soloCBwhitelist 3M-february-2018.txt.gz \
                                --genomeDir ../Index_TAIR11_STAR --readFilesIn $File2 $File1 --outFileNamePrefix ${name}_star \
                                --outSAMtype BAM --soloUMIlen 12 --soloOutFileNames StarCounts_${name}
}

        


# Potential error:
# EXITING because of FATAL ERROR in input read file: the total length of barcode sequence is 26 not equal to expected 28
