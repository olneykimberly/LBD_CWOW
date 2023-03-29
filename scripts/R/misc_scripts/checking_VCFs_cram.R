library(vroom)
setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R")

# read in metadata
pathToRawData = c("/research/labs/neurology/fryer/projects/LBD_CWOW/")
metadata <- vroom(paste0(pathToRawData, "RNA_metadata.tsv"))
WGS_sampleID <- read.delim("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R/CWOW_WGS_230301.txt", sep = "\t", header = T)

shared <- intersect(WGS_sampleID$NPID,metadata$NPID)

WGS_sampleID <- read.delim("/research/labs/neurology/fryer/projects/LBD_CWOW/WGS_sampleID.txt", sep = "\t", header = T)

#----------------
# WGS cram location
# /research/bsi/archive/PI/Ross_Owen_oar01/tertiary/s212492.LBD_project/whole_genome/cram_files_from_NIH
# check cram folder as well. 
cram_files <- list.files("/research/bsi/archive/PI/Ross_Owen_oar01/tertiary/s212492.LBD_project/whole_genome/cram_files_from_NIH/", ".cram.crai")
# remove the .cram.crai from file name to compare to metadata list
cram_files <- gsub(".cram.crai","",cram_files)

# cram file but not in VCF
meta_and_cram <- intersect(metadata$NIH.WGS.ID, cram_files)
inRNAbutNotinCram <- setdiff(metadata$NIH.WGS.ID, cram_files)
inCrambutNotinRNA <- setdiff(cram_files, metadata$NIH.WGS.ID)

# Yes cram file present
WGS_sampleID$cram_present <- ifelse(WGS_sampleID$NIH.WGS.ID %in% cram_files, "yes", "no")
# Yes RNA present
WGS_sampleID$RNA_present <- ifelse(WGS_sampleID$NPID %in% metadata$NPID, "yes", "no")


write.table(WGS_sampleID, "/research/labs/neurology/fryer/projects/LBD_CWOW/WGS_sampleID_cram_present.txt", row.names = F, quote = F, sep = "\t")
#----------------
# AllSamples.variants.vcf.gz
# /research/bsi/archive/PI/Ross_Owen_oar01/tertiary/s212492.LBD_project/whole_genome/joint_vcf/AllSamples.variants.vcf.gz
joint_vcf_sampleIDs <- read.delim("/research/labs/neurology/fryer/projects/LBD_CWOW/joint_vcf_sampleIDs.txt", sep = "\t", header = F)
joint_vcf_sampleIDs <- joint_vcf_sampleIDs[10:961]
joint_vcf_sampleIDs_trans <- as.data.frame(t(joint_vcf_sampleIDs))
names(joint_vcf_sampleIDs_trans)[1] ="NIH.WGS.ID"

# Check joint VCF
meta_and_joint <- intersect(metadata$NIH.WGS.ID, joint_vcf_sampleIDs_trans$NIH.WGS.ID)
inRNAbutNotinjoint <- setdiff(metadata$NIH.WGS.ID, joint_vcf_sampleIDs_trans$NIH.WGS.ID)
injointbutNotinRNA <- setdiff(joint_vcf_sampleIDs_trans$NIH.WGS.ID, metadata$NIH.WGS.ID)

#----------------
# VCFs from NIH
NIH_vcf_sampleIDs <- read.delim("/research/labs/neurology/fryer/projects/LBD_CWOW/vcfs_from_NIH.txt", sep = "\t", header = F)
NIH_vcf_sampleIDs <- NIH_vcf_sampleIDs[10:1796]
NIH_vcf_sampleIDs <- as.data.frame(t(NIH_vcf_sampleIDs))
names(NIH_vcf_sampleIDs)[1] ="NIH.WGS.ID"

meta_and_NIH <- intersect(metadata$NIH.WGS.ID, NIH_vcf_sampleIDs$NIH.WGS.ID)
inRNAbutNotinNIH <- setdiff(metadata$NIH.WGS.ID, NIH_vcf_sampleIDs$NIH.WGS.ID)
inNIHbutNotinRNA <- setdiff(NIH_vcf_sampleIDs$NIH.WGS.ID, metadata$NIH.WGS.ID)

#----------------
# Mayo VCFs
Mayo_vcf_sampleIDs <- read.delim("/research/labs/neurology/fryer/projects/LBD_CWOW/Mayo_VCFs.txt", sep = "\t", header = F)
Mayo_vcf_sampleIDs <- Mayo_vcf_sampleIDs[10:512]
Mayo_vcf_sampleIDs <- as.data.frame(t(Mayo_vcf_sampleIDs))
names(Mayo_vcf_sampleIDs)[1] ="NIH.WGS.ID"

inboth <- intersect(metadata$NIH.WGS.ID, Mayo_vcf_sampleIDs$NIH.WGS.ID)
inRNAbutNotinVCF <- setdiff(metadata$NIH.WGS.ID, Mayo_vcf_sampleIDs$NIH.WGS.ID)
inVCFbutNotinRNA <- setdiff(Mayo_vcf_sampleIDs$NIH.WGS.ID, metadata$NIH.WGS.ID)

