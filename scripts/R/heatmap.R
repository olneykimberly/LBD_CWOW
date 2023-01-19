library(forcats)

ABD$differnce <- ABD$female_logFC - ABD$male_logFC
ABD <- ABD[order(-ABD$differnce),]
TEST <- reshape2::melt(ABD)
data <- subset(TEST, variable == "female_logFC" | variable == "male_logFC")
data$gene_name <- factor(data$gene_name , levels = unique(data$gene_name))
data$gene_name <- fct_rev(data$gene_name)

ggplot(data = data) +
  geom_tile(aes(x = variable, y =gene_name, fill = value)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "rgb",
    guide = "colourbar",
    breaks = c(-1, 0, 1),
    name = expression(log[2](CPM))
  ) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(0, 0.2, 0, 0.2, "cm"), 
    panel.spacing = unit(0,'lines'), 
    plot.title = element_text(size = 10, vjust = -1, hjust = 0),
    axis.title.x=element_blank(),
   # axis.text.x = element_blank(), 
    axis.ticks.x=element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.y = element_text(size = 10, margin = margin(r = -2)),
    axis.ticks.y=element_blank()) +
  ggtitle("DE in XX females") 
