params <-
list(args = "myarg")

#' ---
#' title: "Differential expression in LBD samples"
#' author: "Kimberly Olney"
#' date: "01/20/2022"
#' output:
#'   html_document:
#'     df_print: paged
#'   pdf_document: default
#' params:
#'   args: myarg
#' ---
#' 
#' # Setup
## ----setup----------------------------------------------------------------------------------------------------------------------------------------------------------
# Also do Session > Set Working Directory > Choose Directory
knitr::opts_knit$set(root.dir = ".")

#' 
## ----libraries, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------------
Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
library(rmarkdown)
library(BiocParallel)  # SnowParam()
library(dplyr)  # left_join()
library(edgeR)  # DGEList()
library(limma)  # plotMDS()
library(ggrepel) # geom_text_repel()
library(ggplot2)  # ggplot()
library(gplots)  # heatmap.2()
library(grDevices)  # colorRampPalette()
require(philentropy)  # JSD()
library(rtracklayer)  # import()
library(stringr)  # str_match()
require(variancePartition)  # fitExtractVarPartModel()
library(reshape)  # melt()
library(Glimma)
library(plyr)
library(corrplot)
library("ggpubr")

#' 
#' # User defined variables
## ----set_variables--------------------------------------------------------------------------------------------------------------------------------------------------
control <- "CONTROL"
condition <- "LBD"
control_color <- "gray29"
condition_color <- "red"
myContrasts <- c("LBD - CONTROL")
tool = c("star")
pathToRef = c("/research/labs/neurology/fryer/projects/references/human/")
pathToRawData = c("/research/labs/neurology/fryer/projects/LBD_CWOW/")

#' 
#' # Read data
## ----read_data------------------------------------------------------------------------------------------------------------------------------------------------------
# read in metadata
metadata <- read.delim(paste0(pathToRawData, "metadata.tsv"), 
                       header = TRUE,
                       sep = "\t")

RIN <- read.delim(paste0(pathToRawData, "RIN.tsv"), 
                       header = TRUE,
                       sep = "\t")
RIN$NPID <- RIN$Sample
metadata <- merge(metadata, RIN, by = "NPID")
# read in counts data
counts <- read.delim(
  paste0("../../featureCounts/LBD.counts"),
  header = TRUE,
  sep = "\t", check.names=FALSE)
# get the sampleID from the counts file to match with metadata
sampleIDs <- colnames(counts)
samples <- as.data.frame(sampleIDs)
# split the sampleID by _ 
sample_count_id <- as.data.frame(str_split_fixed(samples$sampleID, "_", 2))
sample_count_id$counts_id <- samples$sampleIDs
# rename columns
names(sample_count_id)[names(sample_count_id) == 'V1'] <- 'NPID'
names(sample_count_id)[names(sample_count_id) == 'V2'] <- 'flow_lane'
# merge sample_count_id and metadata files
# this way metadata contains the sampleID to match the counts table 
counts_metadata <- merge(metadata, sample_count_id, by = "NPID")

# gene information
counts.geneid <- read.delim(
  "../../featureCounts/gene_info.txt",
  header = TRUE,
  sep = "\t"
)
counts.geneid <- as.data.frame(counts.geneid$Geneid)
colnames(counts.geneid) <- "gene_id"
# add gene_name as row names to counts file
rownames(counts) <- counts.geneid$gene_id

# read in annotation file
gtf.file <- paste0(pathToRef, "gencode.v38.annotation.gtf")
genes.gtf <- rtracklayer::import(gtf.file)
genes.gtf <- as.data.frame(genes.gtf)
genes.gtf <- genes.gtf[genes.gtf$type == "gene",]
table(genes.gtf$gene_type)

#counts.geneid
#names(counts.geneid)[1] <- 'gene_name'
joined.df <- join(counts.geneid, genes.gtf, type = "inner")

# check columns and rows match up between files
all.equal(rownames(counts), counts.geneid$gene_id)
all.equal(rownames(counts), genes.gtf$gene_id)

# reorder counts to be in the same order as metadata table
counts <- counts[, counts_metadata$counts_id]
all.equal(colnames(counts), (counts_metadata$counts_id))

#' 
#' # Create DGE object
## ----DGE_object-----------------------------------------------------------------------------------------------------------------------------------------------------
# create object
dge <- DGEList(counts = counts,
               samples = counts_metadata,
               genes = genes.gtf)

table(dge$samples$TYPE)

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# remove AD and PA controls
remove <- vector()
for (i in 1:nrow(dge$samples)) {
  if (dge$samples$TYPE[i] == "CONTROL - AD" | dge$samples$TYPE[i] == "CONTROL - PA") {
    remove <- c(remove, i)
  }
}
dge <- dge[,-remove, keep.lib.sizes = FALSE]
table(dge$samples$TYPE)

#' 
#' # Remove mitochondrial genes
## ----MT_genes-------------------------------------------------------------------------------------------------------------------------------------------------------
dim(dge)
removeMT <- dge$genes$seqnames != "chrM"  # true when NOT MT
dge <- dge[removeMT,,keep.lib.sizes = FALSE]
dim(dge)

#' 
#' # Save functions
#' These functions with help simultaneously save plots as a png, pdf, and tiff 
#' file.
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveToPDF <- function(...) {
    d = dev.copy(pdf,...)
    dev.off(d)
}

saveToPNG <- function(...) {
    d = dev.copy(png,...)
    dev.off(d)
}

#' 
#' # Raw MDS with technical replicates
## ----MDS_techreps, warning=FALSE------------------------------------------------------------------------------------------------------------------------------------
# set colors and get data
table(dge$samples$TYPE)
dge$samples$TYPE <- as.factor(dge$samples$TYPE)
group_colors <- c("grey","purple")[dge$samples$TYPE]
lcpm <- edgeR::cpm(dge$counts, log = TRUE)
par(bg = 'white')

# plot MDS
plotMDS(
  lcpm,
  top = 100, 
  labels = dge$samples$Sex,
  cex = 2, 
  dim.plot = c(1,2), 
  plot = TRUE, 
  col = group_colors, gene.selection = "common"
)
title(expression('Top 100 Genes (Log'[2]~'CPM)'))

path <- paste0("../../results/", tool, "/MDS/", condition,"_gene_MDS_techreps_label_sex")
saveToPDF(paste0(path, ".pdf"), width = 5.2, height = 5.2)
saveToPNG(paste0(path, ".png"), width = 4, height = 4, unit = "in", res = 300)

#' # Sex check
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
genes_and_counts <- cbind(dge$genes$gene_name, lcpm)
genes_and_counts <- as.data.frame(genes_and_counts)
names(genes_and_counts)[names(genes_and_counts) == "V1"] <- "Geneid"
rownames(genes_and_counts)<-NULL
genes_counts <- melt(genes_and_counts, id=c("Geneid"))
names(genes_counts)[names(genes_counts) == "variable"] <- "counts_id"

df <- cbind(counts_metadata$counts_id, counts_metadata$Sex)
df <- as.data.frame(df)
names(df)[names(df) == "V1"] <- "counts_id"
names(df)[names(df) == "V2"] <- "Sex"

data <- merge(genes_counts, df, by = "counts_id")

sexGenes <- c("DDX3X, DDX3Y")
SelectGenes_counts <-
  subset(
    data,
    Geneid %in% c(
      "DDX3X",
      "DDX3Y",
      "ZFX",
      "ZFY",
      "USP9X",
      "USP9Y",
      "KDM6A",
      "UTY",
      "PCDH11X",
      "PCDH11Y",
      "XIST",
      "SRY"
    )
  )
SelectGenes_counts[, "geneComb"] <- NA
SelectGenes_counts[, "group"] <- NA


SelectGenes_counts$geneComb <-
  ifelse(
    SelectGenes_counts$Geneid == "DDX3X",
    "DDX3X:DDX3Y",
    ifelse(
      SelectGenes_counts$Geneid == "DDX3Y",
      "DDX3X:DDX3Y",
      ifelse(
        SelectGenes_counts$Geneid == "ZFX",
        "ZFX:ZFY",
        ifelse(
          SelectGenes_counts$Geneid == "ZFY",
          "ZFX:ZFY",
          ifelse(
            SelectGenes_counts$Geneid == "USP9X",
            "USP9X:USP9Y",
            ifelse(
              SelectGenes_counts$Geneid == "USP9Y",
              "USP9X:USP9Y",
              ifelse(
                SelectGenes_counts$Geneid == "KDM6A",
                "UTX:UTY",
                ifelse(
                  SelectGenes_counts$Geneid == "UTY",
                  "UTX:UTY",
                  ifelse(
                    SelectGenes_counts$Geneid == "PCDH11X",
                    "PCDH11X:PCDH11Y",
                    ifelse(
                      SelectGenes_counts$Geneid == "PCDH11Y",
                      "PCDH11X:PCDH11Y",
                      ifelse(
                        SelectGenes_counts$Geneid == "XIST",
                        "XIST",
                        ifelse(SelectGenes_counts$Geneid == "SRY", "SRY", "NA")
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

SelectGenes_counts$group <-
  ifelse(
    SelectGenes_counts$geneComb == "DDX3X:DDX3Y",
    1,
    ifelse(
      SelectGenes_counts$geneComb == "ZFX:ZFY",
      4,
      ifelse(
        SelectGenes_counts$geneComb == "USP9X:USP9Y",
        3,
        ifelse(
          SelectGenes_counts$geneComb == "UTX:UTY",
          5,
          ifelse(
            SelectGenes_counts$geneComb == "PCDH11X:PCDH11Y",
            2,
            ifelse(
              SelectGenes_counts$geneComb == "XIST",
              6,
              ifelse(SelectGenes_counts$geneComb == "SRY", 7, "NA")
            )
          )
        )
      )
    )
  )
# Plot
data <- SelectGenes_counts
not_male <- subset(data, counts_id == "5510_FCH5MYFDMXY_L1" | counts_id ==  "5510_FCH5MYFDMXY_L2" )
not_male$Sex <- c()
leg_lab <- "reported sex"
cbPaletteJITTER = c("darkorange", "blue")
geneticSEXgenes_plot <- ggplot(data, aes(x = Geneid, y = value)) +
  geom_jitter(aes(color = Sex, shape = Sex),
              width = 0.25,
              size = 2.0) +
  scale_color_manual(leg_lab, values = cbPaletteJITTER) + # Jitter color palette
  scale_shape_manual(leg_lab, values = c(19, 15)) +
  labs(x = "", y = "", title = "") +
  facet_grid(
    . ~ group + geneComb,
    switch = "x",
    # Moves the labels from the top to the bottom
    #labeller = label_both, # Adds the labels to the year and X variables
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) +
  xlab("") + ylab("lcpm") # Removes the month legend

geneticSEXgenes_plot + theme(
  panel.background = element_rect(
    fill = "white",
    colour = "white",
    size = 0.5,
    linetype = "solid"
  ),
  panel.grid.major = element_line(
    size = 0.5,
    linetype = 'solid',
    colour = "white"
  ),
  panel.grid.minor = element_line(
    size = 0.25,
    linetype = 'solid',
    colour = "grey"
  ),
  panel.border = element_rect(
    colour = "black",
    fill = NA,
    size = 1
  ),
  axis.text.x = element_text(face = "italic")
)

path <- paste0("../../results/", tool, "/sex_check/LBD_gene_raw_sex_check")
saveToPDF(paste0(path, ".pdf"), width = 10.5, height = 6)
saveToPNG(paste0(path, ".png"), width = 10.5, height = 4, unit = "in", res = 300)

#' # RIN check with replicates
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# first filter by expression and normalize the data
keep.expr <- filterByExpr(dge, group = dge$samples$TYPE)
dim(dge)
dge.filtered.reps <- dge[keep.expr, , keep.lib.sizes = FALSE]

dim(dge.filtered.reps)
table(dge.filtered.reps$genes$gene_type)

# Now, normalization by the method of trimmed mean of M-values (TMM)
dge.filtered.reps.norm <- calcNormFactors(dge.filtered.reps, method = "TMM")

# norm factor summary
summary(dge.filtered.reps.norm$samples$norm.factors)
log2cpm.norm.reps <- edgeR::cpm(dge.filtered.reps.norm, log = T)
nsamples <- ncol(dge.filtered.reps.norm)
boxplot(log2cpm.norm.reps, 
        main="Filtered normalized lcpm data with replicates", 
        xlab="RIN", 
        ylab=expression('Counts per gene (Log'[2]~'CPM)'),
        axes=FALSE)
axis(2)
axis(1,at=1:nsamples,labels=(dge.filtered.reps.norm$samples$RIN),las=2,cex.axis=0.8)

path <- paste0("../../results/", tool, "/boxplot/LBD_gene_boxplot_with_reps_with_RIN")
saveToPDF(paste0(path, ".pdf"), width = 35, height = 6)
saveToPNG(paste0(path, ".png"), width = 35, height = 6, unit = "in", res = 300)

#' 
#' # Sum technical replicates
## ----techReps-------------------------------------------------------------------------------------------------------------------------------------------------------
# sum technical replicates
dim(dge)
dge.tech <- sumTechReps(dge, dge$samples$NPID)
dim(dge.tech$counts)
colnames(dge.tech$counts) <- dge.tech$samples$NPID

#' 
#' # Raw MDS
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# set colors and get data
group_colors <- c("black", "red")[dge.tech$samples$TYPE]
lcpm <- edgeR::cpm(dge.tech$counts, log = TRUE)

par(bg = 'white')

# plot MDS
plotMDS(
  lcpm, 
  top = 100, 
  labels = dge.tech$samples$NPID,
  cex = .8, 
  dim.plot = c(4,5), 
  plot = TRUE, 
  col = group_colors,
  gene.selection = "common"
)
title(expression('Top 100 Genes - Raw (Log'[2]~'CPM)'))
legend(
  "top",
  legend = c(control, condition),
  pch = 16,
  col = c(control_color, condition_color),
  cex = 1
)
# save
path <- paste0("../../results/", tool, "/MDS/", condition, "_gene_MDS_raw")
saveToPDF(paste0(path, ".pdf"), width = 4, height = 4)
saveToPNG(paste0(path, ".png"), width = 4, height = 4, unit = "in", res = 300)

#' 
#' # Filter lowly expressed genes
#' 
#' The filterByExpr() function in the edgeR package determines which genes have a 
#' great enough count value to keep.  We will filter by group.  This means at least 
#' 6 samples (6 is the smallest group sample size) must express a minimum count of 
#' 10 (in cpm, default value).
#' 
## ----filter---------------------------------------------------------------------------------------------------------------------------------------------------------
keep.expr <- filterByExpr(dge.tech, group = dge.tech$samples$TYPE)
dim(dge.tech)
dge.filtered <- dge.tech[keep.expr, , keep.lib.sizes = FALSE]

dim(dge.filtered)
table(dge.filtered$genes$gene_type)

#' 
#' # TMM normalization
## ----TMM_normalize--------------------------------------------------------------------------------------------------------------------------------------------------
# Now, normalization by the method of trimmed mean of M-values (TMM)
dge.filtered.norm <- calcNormFactors(dge.filtered, method = "TMM")

# norm factor summary
summary(dge.filtered.norm$samples$norm.factors)

#' 
#' # Density plot
#' Density plots of log-intensity distribution of each library can be superposed 
#' on a single graph for a better comparison between libraries and for 
#' identification of libraries with weird distribution. 
## ----density_plots--------------------------------------------------------------------------------------------------------------------------------------------------
# set graphical parameter
par(mfrow = c(1,3))

# Normalize data for library size and expression intesntiy
log2cpm.tech <- edgeR::cpm(dge.tech, log = TRUE)
log2cpm.filtered <- edgeR::cpm(dge.filtered, log = TRUE)
log2cpm.norm <- edgeR::cpm(dge.filtered.norm, log = TRUE)

# set colors
colors <- c("red","orange","green","yellow","blue","purple", 
            "lightgray","brown","pink","cyan")
nsamples <- ncol(dge.tech)

# First, plot the first column of the log2cpm.tech density
plot(density(log2cpm.tech[,1]), col = colors[1], lwd = 2, ylim = c(0,0.25), 
     las = 2, main = "A. Raw", xlab = expression('Log'[2]~CPM))

# For each sample plot the lcpm density
for (i in 2:nsamples){
  den <- density(log2cpm.tech[,i]) #subset each column
  lines(den$x, den$y, col = colors[i], lwd = 2) 
}

# Second, plot log2cpm.filtered
plot(density(log2cpm.filtered[,1]), col = colors[1], lwd = 2, ylim = c(0,0.25), 
     las = 2, main = "B. Filtered", xlab = expression('Log'[2]~CPM))
abline(v = edgeR::cpm(3, log = TRUE), lty = 3)
for (i in 2:nsamples) {
  den <- density(log2cpm.filtered[,i])
  lines(den$x, den$y, col = colors[i], lwd = 2)
}

# Third, plot log2cpm.norm
plot(density(log2cpm.norm[,1]), col = colors[1], lwd = 2, ylim = c(0,0.25), 
     las = 2, main = "C. Normalized", xlab = expression('Log'[2]~CPM))
abline(v = edgeR::cpm(3, log = TRUE), lty = 3)
for (i in 2:nsamples) {
  den <- density(log2cpm.norm[,i])
  lines(den$x, den$y, col = colors[i], lwd = 2)
}

# save
path <- paste0("../../results/", tool, "/density/LBD_gene_density")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)
saveToPNG(paste0(path, ".png"), width = 6, height = 4, unit = "in", res = 300)

#' 
#' # Boxplots
## ----boxplots-------------------------------------------------------------------------------------------------------------------------------------------------------
# set parameters
par(mfrow = c(1,3))

# First look at dge.tech
boxplot(log2cpm.tech, 
        main="A. Raw", 
        xlab="", 
        ylab=expression('Counts per gene (Log'[2]~'CPM)'),
        axes=FALSE,
        col = colors
        )
axis(2) # 2 = left 
axis(1, # 1 = below 
     at = 1:nsamples, # points at which tick-marks should be drawn
     labels = colnames(log2cpm.tech),
     las = 2,
     cex.axis = 0.8 # size of axis
     )

# Second, look at dge.filtered
boxplot(log2cpm.filtered, 
        main="B. Filtered", 
        xlab="", 
        ylab=expression('Counts per gene (Log'[2]~'CPM)'),
        axes=FALSE,
        col = colors
        )
axis(2)
axis(1, at=1:nsamples,labels=colnames(log2cpm.filtered),las=2,cex.axis=0.8)

# Third, look at dge.norm
boxplot(log2cpm.norm, 
        main="C. Normalized", 
        xlab="", 
        ylab=expression('Counts per gene (Log'[2]~'CPM)'),
        axes=FALSE,
        col = colors)
axis(2)
axis(1,at=1:nsamples,labels=colnames(log2cpm.norm),las=2,cex.axis=0.8)

# save
path <- paste0("../../results/", tool, "/boxplot/LBD_gene_boxplot")
saveToPDF(paste0(path, ".pdf"), width = 10, height = 6)
saveToPNG(paste0(path, ".png"), width = 10, height = 6, unit = "in", res = 300)

#' 
#' # RIN check with replicates summed
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
boxplot(log2cpm.norm, 
        main="Filtered normalized lcpm data", 
        xlab="RIN", 
        ylab=expression('Counts per gene (Log'[2]~'CPM)'),
        axes=FALSE,
        col = colors)
axis(2)
axis(1,at=1:nsamples,labels=(dge.filtered.norm$samples$RIN),las=2,cex.axis=0.8)

path <- paste0("../../results/", tool, "/boxplot/LBD_gene_boxplot_sum_reps_with_RIN")
saveToPDF(paste0(path, ".pdf"), width = 12, height = 6)
saveToPNG(paste0(path, ".png"), width = 12, height = 6, unit = "in", res = 300)

# check if there is correlation between RIN and library size
box <- dge.filtered.norm$samples
plot(box$RIN, box$lib.size)

cor(box$RIN, box$lib.size, method = c("pearson", "kendall", "spearman"))
cor.test(box$RIN, box$lib.size, method=c("pearson", "kendall", "spearman"))

# is the data normally distrubited 
ggqqplot(box$lib.size, ylab = "library size")
# wt
ggqqplot(box$RIN, ylab = "RIN")
res <- cor.test(box$lib.size, box$RIN, 
                method = "pearson")
res

ggscatter(box, x = "RIN", y = "lib.size", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "library size", xlab = "RIN value") 

path <- paste0("../../results/", tool, "/library/LBD_corr_RIN_lib_size")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 6)
saveToPNG(paste0(path, ".png"), width = 6, height = 6, unit = "in", res = 300)

#' 
#' # Sex check with summed replicates 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
lcpm <- edgeR::cpm(dge.filtered.norm, log = TRUE)
genes_and_counts <- cbind(dge.filtered.norm$genes$gene_name, lcpm)
genes_and_counts <- as.data.frame(genes_and_counts)
names(genes_and_counts)[names(genes_and_counts) == "V1"] <- "Geneid"
rownames(genes_and_counts)<-NULL
genes_counts <- melt(genes_and_counts, id=c("Geneid"))
names(genes_counts)[names(genes_counts) == "variable"] <- "NPID"

df <- cbind(counts_metadata$NPID, counts_metadata$Sex)
df <- as.data.frame(df)
names(df)[names(df) == "V1"] <- "NPID"
names(df)[names(df) == "V2"] <- "Sex"

data <- merge(genes_counts, df, by = "NPID")

sexGenes <- c("DDX3X, DDX3Y")
SelectGenes_counts <-
  subset(
    data,
    Geneid %in% c(
      "DDX3X",
      "DDX3Y",
      "ZFX",
      "ZFY",
      "USP9X",
      "USP9Y",
      "KDM6A",
      "UTY",
      "PCDH11X",
      "PCDH11Y",
      "XIST",
      "SRY"
    )
  )
SelectGenes_counts[, "geneComb"] <- NA
SelectGenes_counts[, "group"] <- NA


SelectGenes_counts$geneComb <-
  ifelse(
    SelectGenes_counts$Geneid == "DDX3X",
    "DDX3X:DDX3Y",
    ifelse(
      SelectGenes_counts$Geneid == "DDX3Y",
      "DDX3X:DDX3Y",
      ifelse(
        SelectGenes_counts$Geneid == "ZFX",
        "ZFX:ZFY",
        ifelse(
          SelectGenes_counts$Geneid == "ZFY",
          "ZFX:ZFY",
          ifelse(
            SelectGenes_counts$Geneid == "USP9X",
            "USP9X:USP9Y",
            ifelse(
              SelectGenes_counts$Geneid == "USP9Y",
              "USP9X:USP9Y",
              ifelse(
                SelectGenes_counts$Geneid == "KDM6A",
                "UTX:UTY",
                ifelse(
                  SelectGenes_counts$Geneid == "UTY",
                  "UTX:UTY",
                  ifelse(
                    SelectGenes_counts$Geneid == "PCDH11X",
                    "PCDH11X:PCDH11Y",
                    ifelse(
                      SelectGenes_counts$Geneid == "PCDH11Y",
                      "PCDH11X:PCDH11Y",
                      ifelse(
                        SelectGenes_counts$Geneid == "XIST",
                        "XIST",
                        ifelse(SelectGenes_counts$Geneid == "SRY", "SRY", "NA")
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

SelectGenes_counts$group <-
  ifelse(
    SelectGenes_counts$geneComb == "DDX3X:DDX3Y",
    1,
    ifelse(
      SelectGenes_counts$geneComb == "ZFX:ZFY",
      4,
      ifelse(
        SelectGenes_counts$geneComb == "USP9X:USP9Y",
        3,
        ifelse(
          SelectGenes_counts$geneComb == "UTX:UTY",
          5,
          ifelse(
            SelectGenes_counts$geneComb == "PCDH11X:PCDH11Y",
            2,
            ifelse(
              SelectGenes_counts$geneComb == "XIST",
              6,
              ifelse(SelectGenes_counts$geneComb == "SRY", 7, "NA")
            )
          )
        )
      )
    )
  )
# Plot
data <- SelectGenes_counts

leg_lab <- "reported sex"
cbPaletteJITTER = c("darkorange", "blue")
geneticSEXgenes_plot <- ggplot(data, aes(x = Geneid, y = value)) +
  geom_jitter(aes(color = Sex, shape = Sex),
              width = 0.25,
              size = 3.0) +
  scale_color_manual(leg_lab, values = cbPaletteJITTER) + # Jitter color palette
  scale_shape_manual(leg_lab, values = c(19, 15)) +
  labs(x = "", y = "", title = "") +
  facet_grid(
    . ~ group + geneComb,
    switch = "x",
    # Moves the labels from the top to the bottom
    #labeller = label_both, # Adds the labels to the year and X variables
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) +
  xlab("") + ylab("lcpm") # Removes the month legend

geneticSEXgenes_plot + theme(
  panel.background = element_rect(
    fill = "white",
    colour = "white",
    size = 0.5,
    linetype = "solid"
  ),
  panel.grid.major = element_line(
    size = 0.5,
    linetype = 'solid',
    colour = "white"
  ),
  panel.grid.minor = element_line(
    size = 0.25,
    linetype = 'solid',
    colour = "grey"
  ),
  panel.border = element_rect(
    colour = "black",
    fill = NA,
    size = 1
  ),
  axis.text.x = element_text(face = "italic")
)

path <- paste0("../../results/", tool, "/sex_check/LBD_gene_raw_sex_check_sumreps")
saveToPDF(paste0(path, ".pdf"), width = 10.5, height = 6)
saveToPNG(paste0(path, ".png"), width = 10.5, height = 4, unit = "in", res = 300)

#' 
#' # Variance partition
#' CCA Heatmap
## ----CCA_heatmap----------------------------------------------------------------------------------------------------------------------------------------------------
#form <- ~ TYPE + Age + Sex + ClinicalDx +lib.size + LBD.type +CDLB
form <- ~ TYPE + PathDx + AD.subtype + LBD.type + CDLB + Braak.NFT + Thal.amyloid + MF.SP + MF.NFT + MF.LB + Cing.LB + MF.Amyloid + MF.Tau + Cing.Synuclein + CWOW.Category + VaD  + TDP.type + Brain.wt + ClinicalDx + FHx + Duration + Sex + Age + Race + PMI + APOE + MAPT + GRN + TMEM106b + RIN + Total.RNA.ng 
# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
info <- as.data.frame(dge.filtered.norm$samples)
C = canCorPairs( form, info)

# replace NA with zero
C[is.na(C)] = 0
cor.mtest <- function(C, ...) {
    C <- as.matrix(C)
    n <- ncol(C)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(C[, i], C[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(C)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(C)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(C, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=FALSE 
         )

path <- paste0("../../results/", tool ,"/varpart/LBD_gene_CCA")
saveToPDF(paste0(path, ".pdf"), width = 25, height = 25)
saveToPNG(paste0(path, ".png"), width = 23, height = 25, unit = "in", res = 300)

#' 
#' # Save R object
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dge.filtered.norm, file = paste0("../../rObjects/LBD_dge.filtered.norm.rds"))
dge.filtered.norm <- readRDS(paste0("../../rObjects/LBD_dge.filtered.norm.rds"))

#' 
## ----variance_partition,  message=FALSE, warnings=FALSE, eval=FALSE, echo=FALSE-------------------------------------------------------------------------------------
## register(SnowParam(detectCores())) # work in parallel (takes a while to run)
## 
## # geneExpr: matrix of gene expression values
## # info: information/metadata about each sample
## geneExpr <- as.matrix(dge.filtered.norm$counts)
## info <- as.data.frame(dge.filtered.norm$samples)
## info$RIN <- as.integer(info$RIN)
## 
## # age is usually a continuous so model it as a fixed effect "age"
## # group is categorical, so model them as random effects "(1|group)"
## # note the syntax
## form <- ~ (1 | TYPE) +
## (1 | LBD.type) +
## (1 | CWOW.Category) +
## (1 | TDP.type)  +
## (1 | ClinicalDx) +
## (1 | FHx) +
## (1 | Sex) +
## (1 | Race) +
## (1 | TMEM106b) +Braak.NFT +
## Thal.amyloid +
## MF.SP  +
## MF.LB +
## Cing.LB +
## MF.Amyloid +
## MF.Tau +
## Cing.Synuclein +
## VaD +
## TDP.43 +
## Brain.wt +
## Duration +
## Age +
## PMI +
## RIN
## 
## varPart <- fitExtractVarPartModel(geneExpr, form, info)
## vp <- sortCols(varPart)

#' 
## ----variance_violins, message=FALSE, warnings=FALSE, eval=FALSE, echo=FALSE----------------------------------------------------------------------------------------
## # plot
## plotVarPart(vp)
## 
## # save
## path <- paste0("../../results/", tool, "/varpart/LBD_gene_varpart_violins")
## saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)
## saveToPNG(paste0(path, ".png"), width = 6, height = 4, unit = "in", res = 300)

#' 
## ----variance_percent_bars, message=FALSE, warnings=FALSE, eval=FALSE, echo=FALSE-----------------------------------------------------------------------------------
## # plot
## plotPercentBars(vp[1:10,])
## # make a column called gene_id to later merge with genes.gtf
## vp$gene_id <- row.names(vp)
## # dataframe of gene_id and gene_name
## gene_id_name <- cbind(genes.gtf$gene_id, genes.gtf$gene_name)
## test <- as.data.frame(test)
## names(test)[names(test) == 'V1'] <- 'gene_id'
## names(test)[names(test) == 'V2'] <- 'gene_name'
## 
## test2 <- merge(vp, test)
## row.names(test2) <- make.names(test2$gene_name, unique = TRUE)
## test2$gene_id <- NULL
## test2$gene_name <- NULL
## 
## plotPercentBars(test2[1:10,])
## 
## # save
## path <- paste0("../../results/", tool, "/varpart/LPS_",tolower(tissue),"_gene_percent_bars")
## saveToPDF(paste0(path, ".pdf"), width = 6, height = 6)
## saveToPNG(paste0(path, ".png"), width = 6, height = 6, unit = "in", res = 300)

#' 
#' # Top variable genes
## ----variance_group,  message=FALSE, warnings=FALSE, eval=FALSE, echo=FALSE-----------------------------------------------------------------------------------------
## varPart.df <- as.data.frame(test2)
## # sort genes based on variance explained by group
## order.varPart.df <- varPart.df[order(varPart.df$RIN, decreasing = TRUE),]
## head(order.varPart.df["RIN"], 10)
## 
## # sort genes based on variance explained by Sex
## order.varPart.df <- varPart.df[order(varPart.df$Sex, decreasing = TRUE),]
## head(order.varPart.df["Sex"], 10)
## 
## # sort genes based on variance explained by Race
## order.varPart.df <- varPart.df[order(varPart.df$Race, decreasing = TRUE),]
## head(order.varPart.df["Race"], 10)

#' 
#' # Design matrix
## ----design_matrix--------------------------------------------------------------------------------------------------------------------------------------------------
age <- as.numeric(dge.filtered.norm$samples$Age)
RIN <- as.numeric(dge.filtered.norm$samples$RIN)

sex <- as.factor(dge.filtered.norm$samples$Sex)
race <- as.factor(dge.filtered.norm$samples$Race)

group <- interaction(dge.filtered.norm$samples$TYPE)

design <- model.matrix(~ 0 + group + sex + (1|RIN))
colnames(design) <- c(control,condition, "sex", "RIN")

design

#' 
#' # Voom
## ----voom-----------------------------------------------------------------------------------------------------------------------------------------------------------
# voom transform counts
v <- voomWithQualityWeights(dge.filtered.norm,
                            design,
                            plot = TRUE)

# save
path <- paste0("../../results/", tool, "/voom/LBD_gene_mean_var_weights")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)
saveToPNG(paste0(path, ".png"), width = 6, height = 4, unit = "in", res = 300)

# fits linear model for each gene given a series of arrays
fit <- lmFit(v, design)

# contrast design for differential expression
contrasts <- makeContrasts(
  title = myContrasts,  # myContrasts was user input from beginning
  levels = colnames(design))
head(contrasts)

# save contrast names
allComparisons <- colnames(contrasts)
allComparisons # check

# run contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)

# Compute differential expression based on the empirical Bayes moderation of the
# standard errors towards a common value.
veBayesFit <- eBayes(vfit)
plotSA(veBayesFit, main = "Final Model: Mean-variance Trend")

# save
path <- paste0("../../results/", tool, "/voom/LBD_gene_final_mean_var")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)
saveToPNG(paste0(path, ".png"), width = 6, height = 4, unit = "in", res = 300)

#' 
#' # Voom MDS Plot
## ----MDS_voom-------------------------------------------------------------------------------------------------------------------------------------------------------
group_colors <- c("black", "red")[v$targets$TYPE]
names <- v$targets$NPID

plotMDS(
  v, 
  top = 100, 
  labels = names,
  cex = 1, 
  dim.plot = c(1,2), 
  plot = TRUE, 
  col = group_colors
)

title(expression('Top 100 Genes - Voom (Log'[2]~'CPM)'))

# save
# save
path <- paste0("../../results/", tool, "/MDS/LBD_gene_MDS_normalized")
saveToPDF(paste0(path, ".pdf"), width = 5, height = 5)
saveToPNG(paste0(path, ".png"), width = 4, height = 4, unit = "in", res = 300)

#' 
## ----save_EList-----------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(v, file = paste0("../../rObjects/LBD_gene_voom.rds"))
v <- readRDS(paste0("../../rObjects/LBD_gene_voom.rds"))

#' 
#' # Number of DEGs
#' Identify number of differential expressed genes.
## ----decide_tests---------------------------------------------------------------------------------------------------------------------------------------------------
pval <- 0.05

sumTable <- 
  summary(decideTests(
    vfit,  # object
    adjust.method = "BH", # by default the method = "separate"
    p.value = pval,
    lfc = 0  # numeric, minimum absolute log2-fold change required
  ))

print(paste0(" FDRq < ", pval))
sumTable

#' 
#' # Output DEG tables
## ----output_DEG_tables----------------------------------------------------------------------------------------------------------------------------------------------
coef <- 1

for (i in allComparisons) {
  # p < 1, log2fc > 0 
  vTopTableAll <-
    topTable(
      veBayesFit, 
      coef = coef,  
      n = Inf, 
      p.value = 1,
      lfc = 0 
    )
  path <- paste("../../results/", tool, "/DEGs/LBD_gene_DEGs_FDRq1.00.txt", sep = "") 
  write.table(
    vTopTableAll,
    path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # p < 0.05, log2fc > 0
  vTopTable1 <-
    topTable( 
      veBayesFit,  
      coef = coef,  
      n = Inf, 
      p.value = 0.05,
      lfc = 0
    )
  path <- paste("../../results/", tool, "/DEGs/LBD_gene_DEGs_FDRq0.05.txt", sep = "") 
  write.table(
    vTopTable1,
    path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # increment -----------------------------------------------------------------
  coef <- coef + 1
}

#' 
#' Read and save table with all genes (FDRq = 1).
## ----read_DEG_table-------------------------------------------------------------------------------------------------------------------------------------------------
condition_vs_control <- read.table(
  paste0("../../results/", tool, "/DEGs/LBD_gene_DEGs_FDRq1.00.txt", sep = ""),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE)

saveRDS(condition_vs_control, file = paste0("../../rObjects/LBD_gene_table.rds"))

#' 
#' # Assign colors
#' Assign colors  values based on FDRq cutoff of 0.05.
## ----assign_colors--------------------------------------------------------------------------------------------------------------------------------------------------
color_values <- vector()
max <- nrow(condition_vs_control)

for(i in 1:max){
  if (condition_vs_control$adj.P.Val[i] < 0.05){
    if (condition_vs_control$logFC[i] > 0){
      color_values <- c(color_values, 1) # 1 when logFC > 0 and FDRq < 0.05
    }
    else if (condition_vs_control$logFC[i] < 0){
      color_values <- c(color_values, 2) # 2 when logFC < 0 and FDRq < 0.05
    }
  }
  else{
    color_values <- c(color_values, 3) # 3 when FDRq >= 0.05
  }
}

condition_vs_control$color_p0.05 <- factor(color_values)

#' 
#' # Subset genes to label
#' Subset the top 10 up and down-regulated genes
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
up <- condition_vs_control[condition_vs_control$color_p0.05 == 1,]
up10 <- up[1:10,]

down <- condition_vs_control[condition_vs_control$color_p0.05 == 2,]
down <- subset(down, down$logFC < -1.5)
down10 <- down[1:7,]

#' 
#' # Volcano plot
## ----volcano--------------------------------------------------------------------------------------------------------------------------------------------------------
hadjpval <- (-log10(max(
  condition_vs_control$P.Value[condition_vs_control$adj.P.Val < 0.05], 
  na.rm=TRUE)))

p_vol <-
  ggplot(data = condition_vs_control, 
         aes(x = logFC,  # x-axis is logFC
             y = -log10(P.Value),  # y-axis will be -log10 of P.Value
             color = color_p0.05)) +  # color is based on factored color column
  geom_point(alpha = 1.5, size = 1.7) +  # create scatterplot, alpha makes points transparent
  theme_bw() +  # set color theme
  theme(legend.position = "none") +  # no legend
  scale_color_manual(values = c("red", "blue","grey")) +  # set factor colors
  labs(
    title = "", # no main title
    x = expression(log[2](FC)), # x-axis title
    y = expression(-log[10] ~ "(" ~ italic("p") ~ "-value)") # y-axis title
  ) +
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  geom_hline(yintercept = hadjpval,  #  horizontal line
                     colour = "#000000",
                     linetype = "dashed") +
  ggtitle("LBD vs Control\nFDRq < 0.05") +
  theme(plot.title = element_text(size = 10)) +
  geom_text_repel(data = up10,
                  aes(x = logFC, y= -log10(P.Value), label = gene_name), 
                  color = "maroon", 
                  fontface="italic",
                  size = 4, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30)
                  ) +
  geom_text_repel(data = down10,
                  aes(x = logFC, y= -log10(P.Value), label = gene_name), 
                  color = "navyblue", 
                  fontface="italic",
                  size = 4, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15)
                  )  #+
 # scale_y_continuous(breaks = seq(0,8,by=1), limits = c(0,8)) +
 # scale_x_continuous(breaks = seq(-3,3,by=1), limits = c(-3,3))
p_vol

# save
path <- paste0("../../results/", tool, "/volcano/LBD_gene_volcano_FDRq0.05")
saveToPDF(paste0(path, ".pdf"), width = 5.2, height = 5.2)
saveToPNG(paste0(path, ".png"), width = 8, height = 6, unit = "in", res = 300)

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

#' 