# step5_volcano_text
rm(list = ls())
source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")
load("output_step4.RData")
library(ggplot2)
library(ggThemeAssist)

# Volcano plot
x_title <- expression('Log'[2] * ' VHL-OE / Control')
y_title <- expression(paste("Log " * italic("P"), " value"))
x_limit <-
  c(floor(min(dat4_dep$log2FC)), ceiling(max(dat4_dep$log2FC)))
y_limit <-
  c(floor(min(-log10(
    dat4_dep$ad.p.value
  ))), ceiling(max(-log10(
    dat4_dep$ad.p.value
  ))))

volcano <-
  ggplot(data = dat4_dep, aes(
    x = log2FC,
    y = -log10(ad.p.value),
    color = regulation
  )) +
  geom_point(size = 0.5) +
  scale_x_continuous(
    name = x_title,
    breaks = seq(x_limit[1], x_limit[2], 1),
    labels = seq(x_limit[1], x_limit[2], 1),
    limits = x_limit,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = y_title,
    breaks = seq(y_limit[1], y_limit[2], 1),
    labels = seq(y_limit[1], y_limit[2], 1),
    limits = y_limit,
    expand = c(0, 0)
  ) +
  geom_hline(
    yintercept = -log10(p_cut),
    linetype = "dashed",
    color = "grey",
    size = 0.5
  ) +
  geom_vline(
    xintercept = log2(fc_cut),
    linetype = "dashed",
    color = "grey",
    size = 0.5
  ) +
  geom_vline(
    xintercept = log2(1 / fc_cut),
    linetype = "dashed",
    color = "grey",
    size = 0.5
  )

# ggThemeAssistGadget(volcano)
volcano <- volcano + theme(
  plot.subtitle = element_text(vjust = 1),
  plot.caption = element_text(vjust = 1),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black", size = 1),
  axis.title = element_text(colour = "black", face = "bold"),
  axis.text = element_text(colour = "black", face = "bold"),
  plot.title = element_text(colour = "black", face = "bold"),
  panel.background = element_rect(fill = NA),
  legend.position = "none"
)

volcano_2 <- volcano +
  scale_color_manual(values = c("#0000FF", "grey15", "#FF0000"),
                     guide = FALSE)
volcano_2
ggsave(
  filename = "volcano_blue.jpeg",
  plot = volcano_2,
  path = "out_folder",
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 300
)

# label方法,text为要标记基因的accession 号vector
# text <- readClipboard()
if (file.exists("Volcano_label_text.xlsx")) {
  require(openxlsx)
  require(ggrepel)
  text <-
    read.xlsx("Volcano_label_text_2.xlsx",
              sheet = 1,
              colNames = TRUE)
  volcano_text_blue <- volcano_2 +
    geom_text_repel(
      data = subset(dat4_dep, Accession %in% text$Accession),
      aes(label = Symbol),
      size = 2,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  ggsave(
    filename = "volcano_text_blue.jpeg",
    plot = volcano_text_blue,
    path = "out_folder",
    width = 6,
    height = 4.5,
    units = "in",
    dpi = 350
  )
} else{
  stop("we can not find the file:Volcano_label_text.xlsx")
}
# The volcano in GPB didn't show that many gene names, so the code is below
if (file.exists("Volcano_label_text.xlsx")) {
  require(openxlsx)
  require(ggrepel)
  text <-
    read.xlsx("Volcano_label_text.xlsx",
              sheet = 1,
              colNames = TRUE)
  volcano_text_blue <- volcano_2 +
    geom_text_repel(
      data = subset(dat4_dep, Accession %in% text$Accession),
      aes(label = Symbol),
      size = 2,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  ggsave(
    filename = "volcano_text_blue_GPB.jpeg",
    plot = volcano_text_blue,
    path = "out_folder",
    width = 6,
    height = 4.5,
    units = "in",
    dpi = 350
  )
} else{
  stop("we can not find the file:Volcano_label_text.xlsx")
}
