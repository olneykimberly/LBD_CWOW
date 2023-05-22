# load libraries
library(gprofiler2)
library(airway)
library(enrichplot)
library(DOSE) 
library(plyr)
library(scales)
library(forcats)
library(rmarkdown)
library(BiocParallel)
library(dplyr)
library(edgeR)
library(limma)
library(ggrepel)
library(ggplot2)
library(gplots)
library(grDevices)
require(philentropy) #install.packages("philentropy")
library(rtracklayer)
library(stringr)
require(variancePartition) # BiocManager::install("variancePartition")
library(reshape)
library(Glimma)
library(plyr)
library(corrplot)
library(ggpubr)
library(tidyverse)
library(caret)
library(glmnet)
library(vroom)
library(matrixStats)
library("data.table")
library(DESeq2)
library(dittoSeq) #BiocManager::install("dittoSeq")
library(Hmisc)
library(tidyr)
library(gridExtra)
library(grid)
require(openxlsx)
library(UpSetR)
library(mvIC) 
library(RColorBrewer)
# Install information
#library(vctrs, lib.loc = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")
#install.packages("remotes")
#remotes::install_github("GabrielHoffman/mvIC")
#library("llapack")
#remotes::install_github("GabrielHoffman/mvIC")
#BiocManager::install("ShrinkCovMat")
#devtools::install_github("gzt/CholWishart", lib = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")
#install.packages('CholWishart')
#library(wesanderson)
#library(DEGreport)

# To install run:
#devtools::install_github('dviraran/xCell')
# devtools::install_github("hagenaue/BrainInABlender")
library(xCell)
library(devtools)
#library(BrainInABlender) # see https://github.com/hagenaue/BrainInABlender 
source(here::here("/research/labs/neurology/fryer/m239830/tools/BrainInABlender/R/Sir_UnMixALot.R"))
library(reshape2)

# paths, colors, shapes and more
LBD <- "LBD"
AD <- "AD"
PA <- "PA"
CONTROL <- "CONTROL"
control_color <- "#4682B4" # gray
AD_color <- "#B4464B" # yellow gold
PA_color <- "#B4AF46" # brown gold
LBD_color <- "gray35" # green
control_shape <- c(15) # square
AD_shape <- c(16) # circle
PA_shape <- c(17) # triangle
LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B4", "#B4AF46","#B4464B", "gray35")
SexColors <- c("#490092", "#D55E00")
colorbindColors <- dittoColors()
correlationColors <-
  colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))


pathToRef = c("/research/labs/neurology/fryer/projects/references/human/")
pathToRawData = c("/research/labs/neurology/fryer/projects/LBD_CWOW/")

saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}

# read in metadata
# the expanded metadata contains inferred sex, RIN, and WGS sample IDs
# metadata <- read.delim(paste0(pathToRawData, "RNA_metadata.tsv"))
# read in metadata with A-T-S score information
metadata <-
  read.delim(
    "/research/labs/neurology/fryer/m239830/LBD_CWOW/rObjects/metadata_A-T-S_scores.txt"
  )
# keep NPID and ATS information to merge with BinB metadata
#keep <- c("NPID", "A", "T", "S", "ATS", "ATS_names")
#df <- metadata_ATS[, (names(metadata_ATS) %in% keep)]

# Note that BinB cell type information is made in cell_type_enrichment script
#metadata <-
#  read.delim(
 #   "/research/labs/neurology/fryer/m239830/LBD_CWOW/rObjects/metadata_BinB_cellType_Zscore.txt"
 # )

# umbrella cell types rather than primary cell types 
# metadata_BinB_umbrella_cellType_Zscore.txt
metadata <-
  read.delim(
    "/research/labs/neurology/fryer/m239830/LBD_CWOW/rObjects/metadata_BinB_umbrella_cellType_Zscore.txt"
  )
metadata$lib.size.2 <- NULL
metadata$norm.factors.2 <- NULL
metadata$lib.size.1 <- NULL
metadata$norm.factors.1 <- NULL
#metadata$lib.size <- NULL
#metadata$norm.factors <- NULL
# set factor levels
metadata$TYPE <-
  factor(metadata$TYPE, levels = c("CONTROL", "PA", "AD", "LBD"))
metadata$ATS_names <-
  factor(
    metadata$ATS_names,
    levels = c(
      "no pathology",
      "amyloid",
      "amyloid + tau",
      "pure synuclein",
      "low amyloid + synuclein",
      "high amyloid + synuclein",
      "amyloid + synuclein + tau"
    )
  )
