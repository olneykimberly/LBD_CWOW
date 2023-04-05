# Run BIC
setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R")
#-------- set up 
#library('variancePartition')
library('edgeR')
library('BiocParallel')
library('Hmisc')
source(here::here("scripts/R", "file_paths_and_colours.R"))
source(here::here("scripts/R", "gtf_path.R"))
condition <- c("")
tool = c("star")
#-------- read in data 
dge.filtered.norm <- readRDS(paste0("../../rObjects/dge.filtered.norm.rds"))
# some samples are missing RIN values. 
# Replace NA with median RIN. 
# This is necessary to be able include RIN as a covariate in voom
# fill missing values of marks2 with median
dge.filtered.norm$samples$RIN <- impute(dge.filtered.norm$samples$RIN, median)
# one sample is missing VaD information
dge.filtered.norm$samples$VaD <- impute(dge.filtered.norm$samples$VaD, median)
dge.filtered.norm$samples$flowcell_and_lane <- factor(dge.filtered.norm$samples$flowcell_and_lane)
dge.filtered.norm$samples$APOE <- factor(dge.filtered.norm$samples$APOE)

info <- as.data.frame(dge.filtered.norm$samples)
genes <- dge.filtered.norm$genes
log2cpm.norm <- edgeR::cpm(dge.filtered.norm, log = TRUE)
#-------- scale variables
scaled.info <-
  info[c(
    "Race_numeric",
    "RIN",
    "Age",
    "PCT_CODING_BASES",
    "PCT_INTERGENIC_BASES",
    "PCT_INTRONIC_BASES",
    "APOE_E4_allele_count"
  )] %>% scale()
scaled.info.df <- as.data.frame(scaled.info)
# Add scaled information to the metadata called "info"
info_with_scale <- cbind(info, scaled.info.df)

#-------- Voom
formula <- (~ 0 + ATS_names)
voom_with_weights <-
  variancePartition::voomWithDreamWeights(
    counts = dge.filtered.norm$counts,
    formula = formula,
    data = dge.filtered.norm$samples,
    BPPARAM = BiocParallel::SnowParam(cores),
    plot = FALSE
  )
# path <- paste0("../../results/", tool, "/voom/", condition, "_", min_expression, "_gene_mean_voom_with_weights_no_covariates")
#saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)
voomCounts <- voom_with_weights$E

#-------- BIC 
baseFormula <- ~ (1 | ATS_names)
# Combine responses on *rows*
Y = with(
  info,
  rbind(
    sex_inferred,
    Race,
    flowcell_and_lane,
    APOE,
    scaled.info.df$RIN,
    scaled.info.df$Age,
    scaled.info.df$APOE_E4_allele_count,
    scaled.info.df$PCT_CODING_BASES,
    scaled.info.df$PCT_INTERGENIC_BASES,
    scaled.info.df$PCT_INTRONIC_BASES, 
    Astrocyte.Zscore,
    Endothelial.Zscore,
    Microglia.Zscore,
    Mural.Zscore,
    Neuron_All.Zscore,
    Neuron_Interneuron.Zscore, 
    Neuron_Projection.Zscore,
    Oligodendrocyte.Zscore,
    Oligodendrocyte_Immature.Zscore,
    RBC.Zscore,
    Choroid_Plexus.Zscore
  )
)

rownames(Y) <-
  c(
    "sex_inferred",
    "Race",
    "flowcell_and_lane",
    "APOE",
    "RIN",
    "Age",
    "APOE_E4_allele_count",
    "PCT_CODING_BASES",
    "PCT_INTERGENIC_BASES",
    "PCT_INTRONIC_BASES",
    "Astrocyte.Zscore",
    "Endothelial.Zscore",
    "Microglia.Zscore",
    "Mural.Zscore",
    "Neuron_All.Zscore",
    "Neuron_Interneuron.Zscore", 
    "Neuron_Projection.Zscore",
    "Oligodendrocyte.Zscore",
    "Oligodendrocyte_Immature.Zscore",
    "RBC.Zscore",
    "Choroid_Plexus.Zscore"
  )
# variables to consider in the model
# categorical variables must be modeled using (1|)
variables = c(
  "(1|sex_inferred)",
  "(1|Race)",
  "(1|flowcell_and_lane)",
  "(1|APOE)",
  "RIN",
  "Age",
  "APOE_E4_allele_count",
  "PCT_CODING_BASES",
  "PCT_INTERGENIC_BASES",
  "PCT_INTRONIC_BASES",
  "Astrocyte.Zscore",
  "Endothelial.Zscore",
  "Microglia.Zscore",
  "Mural.Zscore",
  "Neuron_All.Zscore",
  "Neuron_Interneuron.Zscore", 
  "Neuron_Projection.Zscore",
  "Oligodendrocyte.Zscore",
  "Oligodendrocyte_Immature.Zscore",
  "RBC.Zscore",
  "Choroid_Plexus.Zscore"
)

# fit forward stepwise regression starting
bestModel_voomcounts = mvForwardStepwise(voomCounts,
                                         baseFormula,
                                         data = info,
                                         variables = variables)
bestModel_voomcounts
