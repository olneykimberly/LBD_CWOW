#!/bin/sh
#$ -cwd
#$ -N r_gene_DE
#$ -m abe
#$ -M olneykimberly@mayo.edu
#$ -l h_vmem=8G
#$ -q 1-day,4-day,lg-mem
#$ -notify


# activate conda environment
source $HOME/.bash_profile
module load python
conda activate LBD

# change directory to where Snakefile is located
CWD="/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R"
cd $CWD

date
Rscript r_gene_DE.r 
date
