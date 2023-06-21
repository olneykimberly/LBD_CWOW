setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/snRNA/scripts/R/")
source(here::here("snRNA/scripts/R", "file_paths_and_colours.R"))
tissue <- "Brain"
treatment <- "pilot"

# read object
dataObject.unannotated <- readRDS("../../rObjects/brain_unannotated.rds")

# set levels and idents
Idents(dataObject.unannotated) <- "seurat_clusters"
# normalize if SCTtransform was used
DefaultAssay(dataObject.unannotated) <- "RNA"
#-------------------------------------------------------------------------------
# Find conserved markers
# initialize df
conserved.markers <- data.frame()
all.clusters <- levels(dataObject.unannotated$SCT_snn_res.0.7)

# loop through each cluster
for (i in all.clusters) {
  
  # print the cluster you're on
  print(i)
  # find conserved marker in cluster vs all other clusters
  markers <- FindConservedMarkers(dataObject.unannotated,
                                  ident.1 = i, # subset to single cluster
                                  grouping.var = "sample", # compare by sample
                                  only.pos = TRUE, # default
                                  min.pct = 0.1, # default
                                  logfc.threshold = 0.25, # default
                                  test.use = "MAST"
  )
  
  # skip if none
  if(nrow(markers) < 3) {
    next
  }
  
  # make rownames a column
  markers <- rownames_to_column(markers, var = "gene")
  
  # make cluster number a column
  markers$cluster <- i
  
  # add to final table
  conserved.markers <- smartbind(conserved.markers, markers)
}
# create delta_pct
conserved.markers$beads_delta_pct <- abs(conserved.markers$beads_pct.1 - 
                                           conserved.markers$beads_pct.2)
conserved.markers$debris_delta_pct <- abs(conserved.markers$debris_pct.1 - 
                                            conserved.markers$debris_pct.2)
# more stringent filtering
markers.strict <- conserved.markers[
  conserved.markers$beads_delta_pct > summary(conserved.markers$beads_delta_pct)[5],]
markers.strict <- conserved.markers[
  conserved.markers$debris_delta_pct > summary(conserved.markers$debris_delta_pct)[5],]
markers.strict$gene_name <- markers.strict$gene
markers.strict$row.num <- 1:nrow(markers.strict)

# compare 
table(as.numeric(conserved.markers$cluster))
table(as.numeric(markers.strict$cluster))

# save (takes a while to run)
saveRDS(conserved.markers, "../../rObjects/brain_conserved_markers.rds")
saveRDS(markers.strict, "../../rObjects/brain_conserved_markers_strict.rds")

#-------------------------------------------------------------------------------
# Find all markers - Finds markers (differentially expressed genes) for each of the identity classes in a dataset
# Find markers for each cluster
# Compares single cluster vs all other clusters
# Default is logfc.threshold = 0.25, min.pct = 0.1
Idents(dataObject.unannotated) <- "seurat_clusters"
all.markers <- FindAllMarkers(object = dataObject.unannotated,
                              assay = "RNA",
                              test.use = "MAST",
                              only.pos = TRUE, # 	Only return positive markers (FALSE by default)
                              verbose = TRUE)
# add column
all.markers$delta_pct <- abs(all.markers$pct.1 - all.markers$pct.2)

# rename columns and rows
rownames(all.markers) <- 1:nrow(all.markers)
all.markers <- all.markers[,c(6,7,1,5,2:4,8)]
colnames(all.markers)[c(6,7)] <- c("pct_1","pct_2")

# save
saveRDS(all.markers, "../../rObjects/unannotated_all_markers.rds")
#-------------------------------------------------------------------------------
# beads vs debris for all clusters before annotation
# initialize variables
clusters <- dataObject.unannotated$seurat_clusters
df <- data.frame()
for (i in clusters) {
  print(i) # inspect 
  # subset object and set idents
  cluster <- subset(dataObject.unannotated, seurat_clusters == i)
  Idents(cluster) <- "sample"
  # differential expression
  markers <- FindMarkers(object = cluster,
                         ident.1 = "beads",
                         ident.2 = "debris",
                         only.pos = FALSE, # default
                         test.use = "MAST",
                         verbose = TRUE,
                         assay = "RNA")
  # if no markers found go to next iteration
  if(nrow(markers) == 0) {
    print(paste0("no markers for ", i))
    next
  }
  # add markers to master df for each comparison
  markers$cluster <- i
  markers$gene <- rownames(markers)
  # combine
  df <- rbind(df, markers)
  # cleanup
  remove(cluster)
}
# reformat table
colnames(df)[c(3,4)] <- c("beads","debris")
rownames(df) <- 1:nrow(df)
df$percent_difference <- abs(df[,3] - df[,4])
df <- df[,c(6,7,1,5,2,3,4,8)]

# write table
path <- "../../results/seurat/DEGs/beads_vs_debris_unannotated.tsv"
write.table(df, path, quote = FALSE, row.names = FALSE, sep = "\t")

