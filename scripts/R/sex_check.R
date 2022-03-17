setwd("~/Dropbox (ASU)/Fryer_Lab/SEPSIS/data_human/")
library(ggplot2)
library(reshape2)
#----------- plot JITTER
#which tissue?

genes_counts <- melt(genes_and_counts)
names(genes_counts)[names(genes_counts) == "variable"] <- "counts_id"

df <- cbind(counts_metadata$counts_id, counts_metadata$Sex)
df <- as.data.frame(df)
names(df)[names(df) == "V1"] <- "counts_id"
names(df)[names(df) == "V2"] <- "Sex"

data <- merge(genes_counts, df, by = "counts_id")
names(data)[names(data) == "genes.gtf$gene_name"] <- "Geneid"

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
data$Geneid <- factor(data$Geneid, as.character(data$Geneid))
data$Geneid
data$geneComb <- factor(data$geneComb, as.character(data$geneComb))
data$Sex <- factor(data$sex, as.character(data$Sex))

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
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) +
  xlab("") # Removes the month legend

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
    colour = "gray"
  ),
  panel.grid.minor = element_line(
    size = 0.25,
    linetype = 'solid',
    colour = "gray"
  ),
  panel.border = element_rect(
    colour = "black",
    fill = NA,
    size = 1
  ),
  axis.text.x = element_text(face = "italic")
)

dev.copy(
  jpeg,
  filename = paste0("figures/sexCheck/", x, "_sexCheckGenes.jpg"),
  width = 10,
  height = 5,
  units = "in",
  res = 1200
)

dev.off()
dev.off()
