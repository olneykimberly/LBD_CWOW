plotSA(fit3)#, xlab="Average log-expression", ylab="log2(sigma)", zero.weights=FALSE, pch=16, cex=0.2, ...)
fit3  <- eBayes(vfit, trend = TRUE)


plotMD(veBayesFit, coef = 1)
veBayesFit$lods

top30 <- order(veBayesFit$lods,decreasing=TRUE)[1:30]
text(veBayesFit$Amean[top30],veBayesFit$coef[top30],labels=veBayesFit$genes[top30,"gene_name"],cex=0.8,col="blue")
topTable(veBayesFit,number=30, coef = 1)


sumTableTest <- decideTests(
    vfit,  # object
    adjust.method = "BH", # by default the method = "separate"
    p.value = pval,
    method = "global",
    lfc = 0  # numeric, minimum absolute log2-fold change required
  )
sumTableTest 
help(decideTests)

test <- (veBayesFit, coef = "PAvsControl")
test <-
  topTable( 
    veBayesFit,  
    coef = "PAvsControl",  
    n = Inf, 
    p.value = 0.05,
    method = "global",
    lfc = 0
  )

type_sex <- interaction(dge.filtered.norm$samples$TYPE, dge.filtered.norm$samples$sex_inferred)
plotMDS(dge.filtered.norm, col = as.numeric(type_sex), gene.selection = "common")

library(factoextra)
library(FactoMineR)

pca.raw.d <- log2(dge.filtered.norm$counts + 0.5)
pca.d <- PCA(t(pca.raw.d), graph = F)
fviz_pca_ind(pca.d, col.ind = type_sex)

library(pheatmap)
design <- model.matrix(~ 0 + TYPE + RIN + Age + sex_inferred + APOE + Race, dge.filtered.norm$samples)
design <- model.matrix(~ 0 + TYPE, dge.filtered.norm$samples)
pheatmap(design)

v <- voom(dge.filtered.norm, design, plot=TRUE)
fit <- lmFit(v, design)
coef.fit <- fit$coefficients

coef = 1
top.table <- topTable(veBayesFit,
n = Inf, 
coef = "ADvsControl",
p.value = 1,
lfc = 0)
head(top.table)
head(top.table[,c(12, 27)])
subset(top.table, gene_name == "AQP1") 
subset(top.table, gene_name == "SERPINA3") 

coef_DEG <- coef.fit[rownames(coef.fit) %in% rownames(top.table)[1:5],]
coef_DEG[,c(1,4)]


