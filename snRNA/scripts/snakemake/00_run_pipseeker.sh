#!/bin/bash
#SBATCH --job-name=pipseeker                         
#SBATCH --partition=cpu-short                          
#SBATCH --nodes=4                                     
#SBATCH --tasks=1                                      
#SBATCH --time=48:00:00 # 8 hours                                
#SBATCH --mem=425G                                        
#SBATCH -o slurm.pipseeker.out
#SBATCH -e slurm.pipseeker.err
#SBATCH --mail-user=olney.kimberly@mayo.edu

# source your bach profile to get your specific settings  
source $HOME/.bash_profile
source $HOME/.bashrc # pipseeker alias is in bashrc 

# change directory
cd /research/labs/neurology/fryer/m239830/LBD_CWOW/snRNA/scripts/snakemake
module load python
conda activate LBD

# 1) get read information
#sh 01_sample_read_info.sh

# 2) create config
#python python 02_create_config.py

# 3) run snakemake - pipseeker alignment via STAR 
snakemake --snakefile 03_pipseeker.Snakefile -j 2 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 32 -t 24:00:00 --mem=120G --partition=cpu-short"

