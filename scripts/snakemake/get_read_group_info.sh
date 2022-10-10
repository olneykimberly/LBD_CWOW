#!/bin/bash

# change directory
cd /research/labs/neurology/fryer/projects/LBD_CWOW/bulkRNA

# create file with list of R1 samples
ls -1 | grep L1_R1_ > R1_samples.txt

# change directory 

# loops through list and print first line
touch sample_read_info.txt
for sample in `cat R1_samples.txt`; do
    #printf "${sample}\t"
    zcat ${sample} | head -1 >> sample_read_info.txt	
done;

mv R1_samples.txt /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/R1_samples.txt 
mv sample_read_info.txt /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/sample_read_info.txt 

cd /research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/snakemake/
paste -d "\t" R1_samples.txt sample_read_info.txt > sample_read_group_info.txt
rm R1_samples.txt
rm sample_read_info.txt

