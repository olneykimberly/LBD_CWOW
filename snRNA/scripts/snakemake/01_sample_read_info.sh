#!/bin/bash

# change directory
cd /research/labs/neurology/fryer/projects/LBD_CWOW/snRNA/

# create file with list of R1 samples
ls -1 | grep _R1 > R1Samples.txt


# loops through list and print first line
touch sampleReadInfo.txt
for sample in `cat R1Samples.txt`; do
    #printf "${sample}\t"
    zcat ${sample} | head -1 >> sampleReadInfo.txt	
done;

mv R1Samples.txt  /research/labs/neurology/fryer/m239830/LBD_CWOW/snRNA/scripts/snakemake/R1Samples.txt
mv sampleReadInfo.txt /research/labs/neurology/fryer/m239830/LBD_CWOW/snRNA/scripts/snakemake/sampleReadInfo.txt

cd /research/labs/neurology/fryer/m239830/LBD_CWOW/snRNA/scripts/snakemake/
paste -d "\t" R1Samples.txt sampleReadInfo.txt > sampleReadGroupInfo.txt
rm R1Samples.txt
rm sampleReadInfo.txt
