# load libraries
library(rmarkdown)
library(BiocParallel)
library(dplyr)
library(edgeR)
library(limma)
library(ggrepel)
library(ggplot2)
library(gplots)
library(grDevices)
require(philentropy)
library(rtracklayer)
library(stringr)
require(variancePartition)
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
library(DEGreport)
library(dittoSeq)
library(Hmisc)
library(mvIC)
library(tidyr)
library(wesanderson)
library(gridExtra)
library(grid)
require(openxlsx)
library(UpSetR)


# To install run:
# devtools::install_github('dviraran/xCell')
# devtools::install_github("hagenaue/BrainInABlender")
library(xCell)
library(devtools)
#library(BrainInABlender) # see https://github.com/hagenaue/BrainInABlender 
source(here::here("scripts/R/brain_blender_functions_modified.R"))
library(reshape2)

# paths, colors, shapes and more
LBD <- "LBD"
AD <- "AD"
PA <- "PA"
CONTROL <- "CONTROL"
control_color <- "#4682B4" # gray
AD_color <- "#B4AF46" # yellow gold
PA_color <- "#B4464B" # brown gold
LBD_color <- "gray35" # green
control_shape <- c(15) # square
AD_shape <- c(16) # circle
PA_shape <- c(17) # triangle
LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B4", "#B4464B", "#B4AF46", "gray35")
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
metadata <- read.delim(paste0(pathToRawData, "RNA_metadata.tsv"))
# read in metadata with A-T-S score information
metadata_ATS <-
  read.delim(
    "/research/labs/neurology/fryer/m239830/LBD_CWOW/rObjects/metadata_A-T-S_scores.txt"
  )
# keep NPID and ATS information to merge with BinB metadata
keep <- c("NPID", "A", "T", "S", "ATS", "ATS_names")
df <- metadata_ATS[, (names(metadata_ATS) %in% keep)]

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