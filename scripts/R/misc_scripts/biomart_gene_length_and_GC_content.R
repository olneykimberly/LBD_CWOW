#BiocManager::install("EDASeq")
library(EDASeq)

# Example:
#ensembl_list <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460")
#getGeneLengthAndGCContent(ensembl_list, "hsa")
#getGeneLengthAndGCContent(ensembl_list, "hsa", mode=c("biomart", "org.db"))

setwd("/Users/m239830/Desktop/")
genes <- read.delim("genes.txt", header = TRUE)
gene_ids <- gsub("\\..*","",genes$gene_id)
genes_1 <- gene_ids[1:500]
genes_2 <- gene_ids[501:1000]
genes_3 <- gene_ids[1001:1500]
genes_4 <- gene_ids[1501:2000]
genes_5 <- gene_ids[2001:2500]
genes_6 <- gene_ids[2500:3000]
genes_7 <- gene_ids[3001:3500]
genes_8 <- gene_ids[3501:4000]
genes_9 <- gene_ids[4001:4500]
genes_10 <- gene_ids[4501:5000]
genes_11 <- gene_ids[5001:5500]
genes_12 <- gene_ids[5501:6000]
genes_13 <- gene_ids[6001:6500]
genes_14 <- gene_ids[6501:7000]
genes_15 <- gene_ids[7001:7500]
genes_16 <- gene_ids[7501:8000]
genes_17 <- gene_ids[8001:8500]
genes_18 <- gene_ids[8501:9000]
genes_19 <- gene_ids[9001:9500]
genes_20 <- gene_ids[9501:10000]
genes_21 <- gene_ids[10001:10500]
genes_22 <- gene_ids[10501:11000]
genes_23 <- gene_ids[11001:11500]
genes_24 <- gene_ids[11501:12000]
genes_25 <- gene_ids[12001:12500]
genes_26 <- gene_ids[12501:13000]
genes_27 <- gene_ids[13001:13500]
genes_28 <- gene_ids[13501:14000]
genes_29 <- gene_ids[14001:14500]
genes_30 <- gene_ids[14501:15000]
genes_31 <- gene_ids[15001:15500]
genes_32 <- gene_ids[15501:16000]
genes_33 <- gene_ids[16001:16500]
genes_34 <- gene_ids[16501:16627]

geneLenght_GCcontent_genes_1 <- getGeneLengthAndGCContent(genes_1, "hsa")
geneLenght_GCcontent_genes_2 <- getGeneLengthAndGCContent(genes_2, "hsa")
geneLenght_GCcontent_genes_3 <- getGeneLengthAndGCContent(genes_3, "hsa")
geneLenght_GCcontent_genes_4 <- getGeneLengthAndGCContent(genes_4, "hsa")
geneLenght_GCcontent_genes_5 <- getGeneLengthAndGCContent(genes_5, "hsa")
geneLenght_GCcontent_genes_6 <- getGeneLengthAndGCContent(genes_6, "hsa")
geneLenght_GCcontent_genes_7 <- getGeneLengthAndGCContent(genes_7, "hsa")
geneLenght_GCcontent_genes_8 <- getGeneLengthAndGCContent(genes_8, "hsa")
geneLenght_GCcontent_genes_9 <- getGeneLengthAndGCContent(genes_9, "hsa")
geneLenght_GCcontent_genes_10 <- getGeneLengthAndGCContent(genes_10, "hsa")

geneLenght_GCcontent_genes_11 <- getGeneLengthAndGCContent(genes_11, "hsa")
geneLenght_GCcontent_genes_12 <- getGeneLengthAndGCContent(genes_12, "hsa")
geneLenght_GCcontent_genes_13 <- getGeneLengthAndGCContent(genes_13, "hsa")
geneLenght_GCcontent_genes_14 <- getGeneLengthAndGCContent(genes_14, "hsa")
geneLenght_GCcontent_genes_15 <- getGeneLengthAndGCContent(genes_15, "hsa")
geneLenght_GCcontent_genes_16 <- getGeneLengthAndGCContent(genes_16, "hsa")
geneLenght_GCcontent_genes_17 <- getGeneLengthAndGCContent(genes_17, "hsa")
geneLenght_GCcontent_genes_18 <- getGeneLengthAndGCContent(genes_18, "hsa")
geneLenght_GCcontent_genes_19 <- getGeneLengthAndGCContent(genes_19, "hsa")
geneLenght_GCcontent_genes_20 <- getGeneLengthAndGCContent(genes_20, "hsa")


geneLenght_GCcontent_genes_21 <- getGeneLengthAndGCContent(genes_21, "hsa")
geneLenght_GCcontent_genes_22 <- getGeneLengthAndGCContent(genes_22, "hsa")
geneLenght_GCcontent_genes_23 <- getGeneLengthAndGCContent(genes_23, "hsa")
geneLenght_GCcontent_genes_24 <- getGeneLengthAndGCContent(genes_24, "hsa")
geneLenght_GCcontent_genes_25 <- getGeneLengthAndGCContent(genes_25, "hsa")
geneLenght_GCcontent_genes_26 <- getGeneLengthAndGCContent(genes_26, "hsa")
geneLenght_GCcontent_genes_27 <- getGeneLengthAndGCContent(genes_27, "hsa")
geneLenght_GCcontent_genes_28 <- getGeneLengthAndGCContent(genes_28, "hsa")
geneLenght_GCcontent_genes_29 <- getGeneLengthAndGCContent(genes_29, "hsa")
geneLenght_GCcontent_genes_30 <- getGeneLengthAndGCContent(genes_30, "hsa")

geneLenght_GCcontent_genes_31 <- getGeneLengthAndGCContent(genes_31, "hsa")
geneLenght_GCcontent_genes_32 <- getGeneLengthAndGCContent(genes_32, "hsa")
geneLenght_GCcontent_genes_33 <- getGeneLengthAndGCContent(genes_33, "hsa")
geneLenght_GCcontent_genes_34 <- getGeneLengthAndGCContent(genes_34, "hsa")

hi <- rbind(geneLenght_GCcontent_genes_1, geneLenght_GCcontent_genes_2, geneLenght_GCcontent_genes_3, geneLenght_GCcontent_genes_4, geneLenght_GCcontent_genes_5, geneLenght_GCcontent_genes_6, geneLenght_GCcontent_genes_7, geneLenght_GCcontent_genes_8, geneLenght_GCcontent_genes_9, geneLenght_GCcontent_genes_10, geneLenght_GCcontent_genes_11, geneLenght_GCcontent_genes_12, geneLenght_GCcontent_genes_13, geneLenght_GCcontent_genes_14, geneLenght_GCcontent_genes_15, geneLenght_GCcontent_genes_17, geneLenght_GCcontent_genes_18, geneLenght_GCcontent_genes_19, geneLenght_GCcontent_genes_20, geneLenght_GCcontent_genes_21, geneLenght_GCcontent_genes_22, geneLenght_GCcontent_genes_23, geneLenght_GCcontent_genes_24, geneLenght_GCcontent_genes_25, geneLenght_GCcontent_genes_26, geneLenght_GCcontent_genes_27, geneLenght_GCcontent_genes_28, geneLenght_GCcontent_genes_29, geneLenght_GCcontent_genes_30, geneLenght_GCcontent_genes_31, geneLenght_GCcontent_genes_32, geneLenght_GCcontent_genes_34, geneLenght_GCcontent_genes_missing)
write.table(hi, "gene_id_gene_length_QC.txt", quote = FALSE, sep = "\t")

# missing
missing_IDs <- rbind(genes_16, genes_33)
write.table(missing_IDs, "missing_IDs.txt", quote = FALSE, sep = "\t")

geneLenght_GCcontent_genes_missing <- getGeneLengthAndGCContent(missing_IDs, "hsa")

