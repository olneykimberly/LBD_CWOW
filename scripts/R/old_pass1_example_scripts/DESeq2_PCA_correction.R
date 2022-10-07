rownames(metadata) <- metadata$RNA.config.ID
help(collapseReplicates)
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~TYPE)

ddsColl <- collapseReplicates(dds, dds$NPID)
dds <- DESeq(dds)
#saveRDS(dds, file = paste0("../../rObjects/LBD_SCC_DESeq.raw.rds"))
summary(dds)
filter_count <- DEGreport::degFilter(counts = counts(dds),
                                     metadata = as.data.frame(colData(dds)),
                                     group = "TYPE",
                                     min = 1, # All samples in group must have more than expr > 0
                                     minreads = 0) 
cat("Genes in final count matrix: ", nrow(filter_count))

vsd <- vst(object = dds, 
           blind = TRUE # blind to experimental design   
)

sampleDists <- vsd %>% 
  assay() %>% 
  t() %>% # Transpose transformed count matrix, such that samples now rows.
  dist(method = "euclidean") # Compute distances between rows of the matrix i.e. sample-to-sample distance

sampleDistMatrix <- sampleDists %>% 
  as.matrix() # Convert to matrix

rownames(sampleDistMatrix) <- str_c(vsd$sample_id, ": ", vsd$TYPE, ", ", vsd$sex_inferred)
colnames(sampleDistMatrix) <- NULL
pheatmap <- pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     col= viridis(n = 20), 
                     main = "Heatmap: uncorrected data")


pca_vsd <- prcomp(t(assay(vsd)))

tibble(PC = c(1:24), 
       sdev = pca_vsd$sdev) %>% 
  dplyr::mutate(d = sdev^2, 
                pev = d/sum(d), 
                cev = cumsum(d)/sum(d))
ggarrange(pcascree(pca_vsd, type = "pev"),
          pcascree(pca_vsd, type = "cev"), 
          nrow = 2,
          labels = c("a", "b"))


list_vsd <- list(pca_vsd)
plot_list <- vector(mode = "list", length = 1)
titles <- c("Uncorrected")

for(i in 1:length(list_vsd)){
  
  plot_list[[i]] <- fviz_pca_ind(list_vsd[[i]],
                                 geom.ind = "point", # show points only (nbut not "text")
                                 col.ind = colData(vsd)[,c("TYPE")], # color by groups
                                 addEllipses = TRUE, ellipse.type = "confidence", # Confidence ellipses
                                 palette =  pal_jco("default", alpha = 0.9)(10)[c(3,9,1,2)],
                                 mean.point = FALSE,
                                 legend.title = "Disease group", 
                                 title = titles[i]
  )
  
}

ggarrange(plotlist = plot_list,
          labels = c("a"),
          nrow = 1)
