#!/bin/bash

for i in *1.fastq.gz; do 
	R1=$i
  # echo $R1 ${R1//_1/_S1_L001_R1_001};
	mv $R1 ${R1//_1/_S1_L001_R1_001};

	#echo ${R1//1.fastq.gz/2.fastq.gz} ${R1//_1/_S1_L001_R2_001};
  mv echo ${R1//1.fastq.gz/2.fastq.gz} ${R1//_1/_S1_L001_R2_001};
	
  # echo ${R1//1.fastq.gz/3.fastq.gz} ${R1//_1/_S1_L001_I1_001};
  mv ${R1//1.fastq.gz/3.fastq.gz} ${R1//_1/_S1_L001_I1_001};

done;
