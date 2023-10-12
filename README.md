# LBD_CWOW
Lewy Body Dementia Center With Out Walls (LBD CWOW) RNAseq processing and analysis.
Anterior cingulate cortex tissue samples from the Mayo Clinic brain bank were collected for 619 individuals. 

| Disease                   | Count   |
| ------------------------- |:-------:|
| Control                   | 86      |
| Pathological  aging (PA)  | 39      |
| Alzheimer’s disease (AD)  | 54      |
| Lewy body dementia (LBD)  | 440     |

This git repo contains scripts for the following:
-   Metadata analysis
-   Processing and analysis bulk-tissue RNA-sequencing data
-   Generation of manuscript figures from Olney et al. 202x publication 
-   Generation of shiny app for exploration of data presented in the the Olney et al. 202x publication, app may be viewed [here](https://fryerlab.shinyapps.io/LBD_CWOW/)


## Create conda environmnet

The necessary software for bulk RNAseq data processing is located here: `LBD.yml`.

To create the environment:
```
conda env create -n LBD --file LBD.yml


# To activate this environment, use
#
#     $ conda activate LBD
#
# To deactivate an active environment, use
#
#     $ conda deactivate

```

#### scripts
snakemake folder contains the scripts for trimming raw fastq files and aligning to a reference genome and generating counts. 
The R folder contains the scripts for differential expression analysis. 

