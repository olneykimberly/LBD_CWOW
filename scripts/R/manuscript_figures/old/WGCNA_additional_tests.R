module.gene.mapping <- as.data.frame(bwnet$colors)
module_colors <- c(unique(bwnet$colors))


module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)
modules_of_interest = c("magenta")

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = voomCounts[genes_of_interest$gene_id,]
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = 12)

row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )
head(edge_list)
library("igraph")
g <- make_empty_graph()
g


plotEigengeneNetworks(bwnet$MEs, setLabels = "purple")
plotNetworkHeatmap(bwnet$MEs, plotGenes = c(bwnet$colors))

purple_graph <- igraph::graph.adjacency(TOM)
plot(purple_graph)
# https://www.biostars.org/p/426343/
#https://www.biostars.org/p/402720/#402733
adj <- TOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network)  
V(network)$color <- bwnet$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
plot(network, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)


#install.packages("remotes")
#remotes::install_github("jtlovell/limmaDE2")
library(limmaDE2)
wgcna2igraph(bwnet$)
igraph::graph.adjacency(bwnet)

counts <- t(voomCounts)
counts[] <- sapply(counts, as.numeric)
soft_power <- 12
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM 
# calculated during module detection, but let us do it again here. 
dissTOM = 1-TOMsimilarityFromExpr(counts, power = 12); 
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap 
plotTOM = dissTOM^7; 
# Set diagonal to NA for a nicer plot diag(plotTOM) = NA; 
# Call the plot function
sizeGrWindow(9,9) 
TOMplot(plotTOM, geneTree, module_colors, main = "Network heatmap plot, all genes")


# https://support.bioconductor.org/p/84412/
g <- graph.adjacency(as.matrix(as.dist(cor(t(MyData), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)

# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)

# Basic chart
plot(network)

# https://r-graph-gallery.com/250-correlation-network-with-igraph.html
# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)

# Basic chart
plot(network)

# https://groups.google.com/g/cytoscape-discuss/c/i7RTJSr7ur4
p1<-adjacency(counts,type=input$adj)

edge<-as.data.frame(get.edgelist(simplify(graph_from_adjacency_matrix(p1, weighted=T))))

#------------------------------------------------------------------------------------
if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install("BioNERO")
library(BioNERO)
set.seed(28) 
setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R")
source(here::here("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R", "file_paths_and_colours.R"))
condition <- c("TYPE") 
tool <- c("star")
dge.filtered.norm <- readRDS(paste0("../../rObjects/dge.filtered.norm.rds"))
info <- as.data.frame(dge.filtered.norm$samples)
df <- subset(info, select = c("TYPE"))
genes <- dge.filtered.norm$genes
voomCounts <- readRDS(paste0("../../rObjects/TYPE.voomCountsMatrix.rds"))
#--------------------
rowRanges <- GRanges(genes)
colData <- DataFrame(info)

final_exp <- SummarizedExperiment(assays=list(counts=voomCounts),
                     rowRanges=rowRanges, colData=df)
help(plot_heatmap)
df = subset(info, select = c("TYPE"))
# Heatmap of sample correlations
p <- plot_heatmap(voomCounts, col_metadata = df,  type = "samplecor", show_rownames = FALSE)
p

# Heatmap of gene expression (here, only the first 50 genes)
p <- plot_heatmap(
  voomCounts[1:50, ], col_metadata = df, type = "expr", show_rownames = FALSE, show_colnames = FALSE
)
p
plot_PCA(final_exp)


net <- exp2gcn(
  final_exp, net_type = "signed hybrid", SFTpower = 12, 
  cor_method = "pearson"
)
# Dendro and colors
plot_dendro_and_colors(net)
# Eigengene networks
plot_eigengene_network(net)
plot_ngenes_per_module(net)


plot_expression_profile(
  exp = final_exp, 
  net = net, 
  plot_module = TRUE, 
  modulename = "yellow"
)