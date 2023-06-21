# load libraries
library("data.table")
library(BiocParallel)
library(DESeq2)
library(DOSE) 
library(DoubletFinder)  
library(Glimma)
library(Hmisc)
library(RColorBrewer)
library(Seurat) 
library(Seurat.utils)
library(SeuratDisk)  
library(ShinyCell)
library(UpSetR)
library(airway)
library(caret)
library(corrplot)
library(cowplot) 
library(devtools)
library(dittoSeq) 
library(dittoSeq) #BiocManager::install("dittoSeq")
library(dplyr)
library(dplyr)  
library(edgeR)
library(enrichplot)
library(forcats)
library(ggplot2)
library(ggplot2) 
library(ggpubr)
library(ggrepel)
library(ggrepel) 
library(ggtree)
library(glmnet)
library(gplots)
library(gprofiler2)
library(grDevices)
library(grid)
library(grid) 
library(gridExtra)
library(gridExtra)  
library(harmony) 
library(limma)
library(matrixStats)
library(mvIC) 
library(parallel)
library(plotly)  
library(plyr)
library(reshape)
library(reshape2)
library(rmarkdown)
library(rtracklayer)
library(scales)
library(shiny)
library(stringr)
library(stringr) 
library(tibble)  
library(tidyr)
library(tidyverse, lib.loc = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")
library(variancePartition)
library(vroom)
require(openxlsx)
require(philentropy) 
require(variancePartition) 
library(gtools)
library(SingleCellExperiment)
library(Matrix.utils)


# paths, colors, shapes and more
LBD <- "LBD"
AD <- "AD"
PA <- "PA"
CONTROL <- "CONTROL"
control_color <- "#4682B4"
AD_color <- "#B4464B" 
PA_color <- "#B4AF46" 
LBD_color <- "gray35" 
control_shape <- c(15) # square
AD_shape <- c(16) # circle
PA_shape <- c(17) # triangle
LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B4", "#B4AF46","#B4464B", "gray35")
SexColors <- c("#490092", "#D55E00")
colorbindColors <- dittoColors()
colors <- colorbindColors
correlationColors <-
  colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

sample_colors <- c("#4682B4", "#B4AF46")
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
