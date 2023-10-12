# LBD_CWOW
Lewy Body Dementia Center With Out Walls (LBD CWOW) RNAseq processing and analysis.
Anterior cingulate cortex tissue samples from the Mayo Clinic brain bank were collected for 619 individuals. 
Raw fastq files are available on SRA, PRJNA1023207

| Disease                   | Count   |
| ------------------------- |:-------:|
| Control                   | 86      |
| Pathological  aging (PA)  | 39      |
| Alzheimer’s disease (AD)  | 54      |
| Lewy body dementia (LBD)  | 440     |

This git repo contains scripts for the following:
-   Metadata analysis
-   Processing of bulk-tissue RNA-sequencing data
-   Analysis of bulk-tissue RNA-sequencing data 
-   Generation of manuscript figures from Olney et al. 202x publication 
-   Generation of shiny app for exploration of the results presented in Olney et al. 202x publication, view app [here](https://fryerlab.shinyapps.io/LBD_CWOW/)


## 1. Create conda enviroment

The necessary software for bulk RNAseq data processing is contained in this file: `LBD.yml`.

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
After the conda LBD environment has been created, you will need to add gatk..
This step must be done manually and not through conda; see [here](https://gatk.broadinstiitute.org
/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4)
The above link will explain how to download GATK4, then you will need to add an alias to your bash profile:
```
alias gatk='/path/to/gatk-package/gatk'
```
Follow the step [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)

