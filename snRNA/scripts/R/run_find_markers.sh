#!/bin/bash
#SBATCH --job-name=FindMarkers                         
#SBATCH --partition=cpu-short                          
#SBATCH --nodes=8                                     
#SBATCH --tasks=1                                      
#SBATCH --time=48:00:00                              
#SBATCH --mem=175G                                        
#SBATCH -o slurm.LBD_markers.out
#SBATCH -e slurm.LBD_markers.err
#SBATCH --mail-user=olney.kimberly@mayo.edu


# activate conda environment
source $HOME/.bash_profile
source ~/.bash_mayobiomods

date
# change directory
cd /research/labs/neurology/fryer/m239830/LBD_COW/snRNA/scripts/R/

#Rscript  BIC_with_cellType_zscores_by_ATS.R 
Rscript find_markers.R

date
