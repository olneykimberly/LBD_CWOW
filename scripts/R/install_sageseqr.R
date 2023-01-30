if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicTools", force = TRUE)

remotes::install_github("DiseaseNeuroGenomics/variancePartition")
install.packages("GenomicTools")

library(snpStats)
library("devtools")
install_github("fischuu/GenomicTools")
library("GenomicTools")
install.packages("gMWT")
remotes::install_github("cran/gMWT")
library(gMWT)
remotes::install_github("fischuu/GenomicTools")
library(GenomicTools)
remotes::install_github("Sage-Bionetworks/sageseqr")
library(sageseqr)
