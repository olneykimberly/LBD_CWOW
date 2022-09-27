setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R/")
dim(dge.filtered.norm)
#removePA <- dge.filtered.norm$samples$TYPE != "CONTROL_PA"  # true when NOT MT
dge.filtered.norm_exPA <- dge.filtered.norm[,which(!dge.filtered.norm$samples$TYPE == "CONTROL_PA")]
dge.filtered.norm_exPA_AD <- dge.filtered.norm_exPA[,which(!dge.filtered.norm_exPA$samples$TYPE == "CONTROL_AD")]

dim(dge.filtered.norm_exPA_AD)
dge.filtered.norm_exPA$samples$TYPE

group_colors <- c("grey50", 
                  control_AD_color, 
                  control_PA_color,
                  LBD_color)[dge.filtered.norm_exPA_AD$samples$TYPE]
Sex <- dge.filtered.norm_exPA_AD$samples$Sex

plotMDS(
  dge.filtered.norm_exPA_AD, 
  top = 100, 
  labels = Sex,
  cex = 1.2, 
  dim.plot = c(1,2), 
  plot = TRUE, 
  col = group_colors, 
  gene.selection = "common"
)
title(expression('Top 100 Genes - Voom (Log'[2]~'CPM)'))
path <- paste0("../../results/", tool, "/MDS/LBD_bothSexes")
saveToPDF(paste0(path, ".pdf"), width = 5, height = 5)

#---- Female only 
dge.filtered.norm_exPA_AD_FEMALE <- dge.filtered.norm_exPA_AD[,which(!dge.filtered.norm_exPA_AD$samples$Sex == "M")]
dim(dge.filtered.norm_exPA_AD_FEMALE)
dge.filtered.norm_exPA_AD_FEMALE$samples$TYPE

group_colors <- c("#D95F02", 
                  control_AD_color, 
                  control_PA_color,
                  LBD_color)[dge.filtered.norm_exPA_AD_FEMALE$samples$TYPE]
Sex <- dge.filtered.norm_exPA_AD_FEMALE$samples$Sex
TYPE <- dge.filtered.norm_exPA_AD_FEMALE$samples$TYPE

plotMDS(
  dge.filtered.norm_exPA_AD_FEMALE, 
  top = 100, 
  labels = TYPE,
  cex = 1.25, 
  dim.plot = c(1,2), 
  plot = TRUE, 
  col = group_colors, 
  gene.selection = "common"
)
title(expression('Top 100 Genes - Voom (Log'[2]~'CPM)'))
path <- paste0("../../results/", tool, "/MDS/LBD_females")
saveToPDF(paste0(path, ".pdf"), width = 5, height = 5)

group <- as.factor(dge.filtered.norm_exPA_AD_FEMALE$samples$TYPE)
design <- model.matrix(~ 0 + group)
colnames(design) <- c("CONTROL", "AD", "PA", "LBD")
v <- voom(dge.filtered.norm_exPA_AD_FEMALE, design, plot=TRUE)
fit <- lmFit(v, design)
contrasts <- makeContrasts(
  LBDvsControl = LBD - CONTROL,
  levels = colnames(design))
head(contrasts)
allComparisons <- colnames(contrasts)
allComparisons # check
vfit <- contrasts.fit(fit, contrasts = contrasts)
veBayesFit <- eBayes(vfit)
plotSA(veBayesFit, main = "Final Model: Mean-variance Trend")
sumTable <- 
  summary(decideTests(
    vfit,  # object
    adjust.method = "BH", # by default the method = "separate"
    p.value = 0.01,
    lfc = 1  # numeric, minimum absolute log2-fold change required
  ))
print(paste0(" FDRq < ", pval))
sumTable
topTable(veBayesFit, coef = "LBDvsControl")

color_values <- vector()
vTopTableAll <-
  topTable(
    veBayesFit, 
    coef = 1,  
    n = Inf, 
    p.value = 1,
    lfc = 0 
  )
female_LBD <- vTopTableAll
max <- nrow(vTopTableAll)
for(i in 1:max){
  if (vTopTableAll$adj.P.Val[i] < 0.01){
    if (vTopTableAll$logFC[i] > 0){
      color_values <- c(color_values, 1) # 1 when logFC > 0 and FDRq < 0.05
    }
    else if (vTopTableAll$logFC[i] < 0){
      color_values <- c(color_values, 2) # 2 when logFC < 0 and FDRq < 0.05
    }
  }
  else{
    color_values <- c(color_values, 3) # 3 when FDRq >= 0.05
  }
}
vTopTableAll$color_p0.05 <- factor(color_values)
up <- vTopTableAll[vTopTableAll$color_p0.05 == 1,]
up10 <- up[1:20,]

down <- vTopTableAll[vTopTableAll$color_p0.05 == 2,]
down10 <- down[1:20,]
hadjpval <- (-log10(max(
  vTopTableAll$P.Value[vTopTableAll$adj.P.Val < 0.01], 
  na.rm=TRUE)))

p_vol <-
  ggplot(data = vTopTableAll, 
         aes(x = logFC,  # x-axis is logFC
             y = -log10(P.Value),  # y-axis will be -log10 of P.Value
             color = color_p0.05)) +  # color is based on factored color column
  geom_point(alpha = 1.5, size = 1.2) +  # create scatterplot, alpha makes points transparent
  theme_bw() +  # set color theme
  theme(legend.position = "none") +  # no legend
  scale_color_manual(values = c("red", "blue", "grey")) +  # set factor colors
  labs(
    title = "", # no main title
    x = expression(log[2](FC)), # x-axis title
    y = expression(-log[10] ~ "(" ~ italic("p") ~ "-value)") # y-axis title
  ) +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  geom_hline(yintercept = hadjpval,  #  horizontal line
             colour = "#000000",
             linetype = "dashed") +
  ggtitle("female (46, XX) LBD vs Control\nFDRq < 0.01") +
  theme(plot.title = element_text(size = 18)) +
  geom_text_repel(data = up10,
                  aes(x = logFC, y= -log10(P.Value), label = gene_name), 
                  color = "maroon", 
                  fontface="italic",
                  size = 5, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15)
  ) +
  geom_text_repel(data = down10,
                  aes(x = logFC, y= -log10(P.Value), label = gene_name), 
                  color = "navyblue", 
                  fontface="italic",
                  size = 5, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10)
  ) +
  scale_y_continuous(breaks = seq(0,15,by=1), limits = c(0,15)) +
  scale_x_continuous(breaks = seq(-2,2.5,by=0.5), limits = c(-2,2.5))
p_vol
path <- paste0("../../results/", tool, "/volcano/female_LBD_gene_volcano_FDRq0.05")
saveToPDF(paste0(path, ".pdf"), width = 7, height = 7)

#---- Male only 
dge.filtered.norm_exPA_AD_MALE <- dge.filtered.norm_exPA_AD[,which(!dge.filtered.norm_exPA_AD$samples$Sex == "F")]
dim(dge.filtered.norm_exPA_AD_MALE)
dge.filtered.norm_exPA_AD_MALE$samples$TYPE

group_colors <- c("#D95F02", 
                  control_AD_color, 
                  control_PA_color,
                  LBD_color)[dge.filtered.norm_exPA_AD_MALE$samples$TYPE]
Sex <- dge.filtered.norm_exPA_AD_MALE$samples$Sex
TYPE <- dge.filtered.norm_exPA_AD_MALE$samples$TYPE

plotMDS(
  dge.filtered.norm_exPA_AD_MALE, 
  top = 100, 
  labels = TYPE,
  cex = 1.25, 
  dim.plot = c(1,2), 
  plot = TRUE, 
  col = group_colors, 
  gene.selection = "common"
)
title(expression('Top 100 Genes - Voom (Log'[2]~'CPM)'))
path <- paste0("../../results/", tool, "/MDS/LBD_males")
saveToPDF(paste0(path, ".pdf"), width = 5, height = 5)

group <- as.factor(dge.filtered.norm_exPA_AD_MALE$samples$TYPE)
design <- model.matrix(~ 0 + group)
colnames(design) <- c("CONTROL", "AD", "PA", "LBD")
v <- voom(dge.filtered.norm_exPA_AD_MALE, design, plot=TRUE)
fit <- lmFit(v, design)
contrasts <- makeContrasts(
  LBDvsControl = LBD - CONTROL,
  levels = colnames(design))
head(contrasts)
allComparisons <- colnames(contrasts)
allComparisons # check
vfit <- contrasts.fit(fit, contrasts = contrasts)
veBayesFit <- eBayes(vfit)
plotSA(veBayesFit, main = "Final Model: Mean-variance Trend")
sumTable <- 
  summary(decideTests(
    vfit,  # object
    adjust.method = "BH", # by default the method = "separate"
    p.value = 0.01,
    lfc = 1  # numeric, minimum absolute log2-fold change required
  ))
print(paste0(" FDRq < ", pval))
sumTable
topTable(veBayesFit, coef = "LBDvsControl")

color_values <- vector()
vTopTableAll <-
  topTable(
    veBayesFit, 
    coef = 1,  
    n = Inf, 
    p.value = 1,
    lfc = 0 
  )
male_LBD <- vTopTableAll
max <- nrow(vTopTableAll)
for(i in 1:max){
  if (vTopTableAll$adj.P.Val[i] < 0.01){
    if (vTopTableAll$logFC[i] > 0){
      color_values <- c(color_values, 1) # 1 when logFC > 0 and FDRq < 0.05
    }
    else if (vTopTableAll$logFC[i] < 0){
      color_values <- c(color_values, 2) # 2 when logFC < 0 and FDRq < 0.05
    }
  }
  else{
    color_values <- c(color_values, 3) # 3 when FDRq >= 0.05
  }
}
vTopTableAll$color_p0.05 <- factor(color_values)
up <- vTopTableAll[vTopTableAll$color_p0.05 == 1,]
up10 <- up[1:20,]

down <- vTopTableAll[vTopTableAll$color_p0.05 == 2,]
down10 <- down[1:20,]
hadjpval <- (-log10(max(
  vTopTableAll$P.Value[vTopTableAll$adj.P.Val < 0.01], 
  na.rm=TRUE)))

p_vol <-
  ggplot(data = vTopTableAll, 
         aes(x = logFC,  # x-axis is logFC
             y = -log10(P.Value),  # y-axis will be -log10 of P.Value
             color = color_p0.05)) +  # color is based on factored color column
  geom_point(alpha = 1.5, size = 1.2) +  # create scatterplot, alpha makes points transparent
  theme_bw() +  # set color theme
  theme(legend.position = "none") +  # no legend
  scale_color_manual(values = c("red", "blue", "grey")) +  # set factor colors
  labs(
    title = "", # no main title
    x = expression(log[2](FC)), # x-axis title
    y = expression(-log[10] ~ "(" ~ italic("p") ~ "-value)") # y-axis title
  ) +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14)) +
  geom_hline(yintercept = hadjpval,  #  horizontal line
             colour = "#000000",
             linetype = "dashed") +
  ggtitle("male (46, XY) LBD vs Control\nFDRq < 0.01") +
  theme(plot.title = element_text(size = 18)) +
  geom_text_repel(data = up10,
                  aes(x = logFC, y= -log10(P.Value), label = gene_name), 
                  color = "maroon", 
                  fontface="italic",
                  size = 5, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15)
  ) +
  geom_text_repel(data = down10,
                  aes(x = logFC, y= -log10(P.Value), label = gene_name), 
                  color = "navyblue", 
                  fontface="italic",
                  size = 5, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15)
  )  +
  scale_y_continuous(breaks = seq(0,15,by=1), limits = c(0,15)) +
  scale_x_continuous(breaks = seq(-1.5,1.5,by=0.5), limits = c(-1.5,1.5))
p_vol
path <- paste0("../../results/", tool, "/volcano/male_LBD_gene_volcano_FDRq0.05")
saveToPDF(paste0(path, ".pdf"), width = 7, height = 7)

write.table(male_LBD, paste0("../../results/", tool, "/DEGs/male_LBD_gene.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
write.table(female_LBD, paste0("../../results/", tool, "/DEGs/female_LBD_gene.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
