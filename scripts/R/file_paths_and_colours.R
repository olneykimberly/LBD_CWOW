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

TypeColors <- c("#4682B4", "#B4464B","#B4AF46", "gray35")
SexColors <- c("orange", "blue")
colorbindColors <- dittoColors()
correlationColors <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))


pathToRef = c("/research/labs/neurology/fryer/projects/references/human/")
pathToRawData = c("/research/labs/neurology/fryer/projects/LBD_CWOW/")

saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

# read in metadata
# the expanded metadata contains inferred sex, RIN, and WGS sample IDs
metadata <- vroom(paste0(pathToRawData, "RNA_metadata.tsv"))
# remove duplicates if any 
metadata <- metadata[!duplicated(metadata[,c('RNA.config.ID')]),]
# rename columns for easier plotting 
# names(metadata)[names(metadata) == "TYPE"] <- "Disease group"
metadata <- metadata %>%
  mutate(across('TYPE', str_replace, 'CONTROL - AD', 'AD'))
metadata <- metadata %>%
  mutate(across('TYPE', str_replace, 'CONTROL - PA', 'PA'))
metadata$TYPE <- factor(metadata$TYPE, levels = c("CONTROL", "AD", "PA", "LBD"))


metadata <- metadata %>% mutate(sex_chr = case_when(
  startsWith(sex_inferred, "female") ~ "XX", 
  startsWith(sex_inferred, "male") ~ "XY"))

# remove rows with no inferred sex.
# No counts data for this samples because there is no RNAseq data for that individual
metadata <- metadata[!is.na(metadata$sex_inferred),]
