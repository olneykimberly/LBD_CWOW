ggplot(data = ABD, aes(x = female_logFC, y = male_logFC, text = paste(gene_name))) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 1,
    ymin = 0.5,
    ymax = 0,
    fill = "lightpink3",
    alpha = .5
  ) +  annotate(
    "rect",
    xmin = 0,
    xmax = -1,
    ymin = 0,
    ymax = -.5,
    fill = "cadetblue3",
    alpha = .5) +
  geom_abline(color = "gray40") +
  geom_text_repel(
    data = LBD_FandM_goi, 
    aes(
      x = female_logFC, 
      y = male_logFC,
      label = gene_name
    ),
    color = "black",
    fontface = "italic",
    size = 2.5, 
    max.overlaps = getOption("ggrepel.max.overlaps", default = 15)
  ) +
  geom_point(size = 1) + 
  theme_bw() +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(size = 10)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  labs(
    title = "DEGs shared and unique between the sexes",
    x = expression(paste("XX female ", log[2](LBD/control))),
    y = expression(paste("XY male ", log[2](LBD/control))))+
  scale_y_continuous(breaks = seq(-1, 1, by = .5), 
                     limits = c(-1, 1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-1, 1, by = .5), 
                     limits = c(-1, 1), expand = c(0,0)) + 
  stat_cor(method = "pearson", label.x = -1.25, label.y = 1)
