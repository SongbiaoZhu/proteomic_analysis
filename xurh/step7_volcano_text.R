# step7_volcano_text
rm(list = ls())
options(stringsAsFactors = F)
source("udf.R")
creDir("data")
creDir("res")
creDir("public")
dataDir = file.path(getwd(),'data')
resDir = file.path(getwd(),'res')
publicDir = file.path(getwd(),'public')
load("output_step1.RData")
load("output_step5.RData")
library(ggplot2)
theme_set(
  theme_gray() + theme(
    axis.line = element_line(size = 0.5),
    panel.background = element_rect(fill = NA, size = rel(20)),
    panel.grid.minor = element_line(colour = NA),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )
)
# Volcano plot
y_title = expression(paste("-Log ", italic("P"), " value"))
volcano <-
  ggplot(data = dat4_dep, aes(
    x = log2FC,
    y = -log10(p.value),
    color = regulation )) + 
  geom_point(size = 0.5) + 
  scale_x_continuous(limits = c(floor(min(dat4_dep$log2FC)/0.5)*0.5,ceiling(max(dat4_dep$log2FC)/0.5)*0.5 ), breaks = seq(floor(min(dat4_dep$log2FC)/0.5)*0.5,ceiling(max(dat4_dep$log2FC)/0.5)*0.5 , 0.5)) +
  scale_y_continuous(limits = c(0, ceiling(max(-log10(dat4_dep$p.value)/0.5)*0.5)), breaks = seq(0,  ceiling(max(-log10(dat4_dep$p.value)/0.5)*0.5),0.5)) +
  xlab("Log2 fold change") + 
  ylab(y_title) +
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
    xintercept = log2(1/fc_cut),
    linetype = "dashed",
    color = "grey",
    size = 0.5
  ) 
# volcano_1 <- volcano + 
#   scale_color_manual(values = c("green", "grey15", "red"),guide=FALSE)
# ggsave(file.path(resDir, "volcano_green.jpeg"),  dpi = 300)
volcano_2 <- volcano + 
  scale_color_manual(values = c("#0000FF", "grey15", "#FF0000"),guide=FALSE)
ggsave(file.path(resDir, "volcano_blue.jpeg"),  dpi = 300)

# label方法,text为要标记基因的accession 号vector
# text <- readClipboard()
if(file.exists(file.path(dataDir, "Volcano_label_text.xlsx"))){
  require(openxlsx)
  require(ggrepel)
  text <- read.xlsx("Volcano_label_text.xlsx", sheet=1, colNames = TRUE )
  # volcano_text_green <- volcano_1 +
  #   geom_text_repel(
  #     data = subset(dat4_dep, Accession %in% text$Accession ),
  #     aes(label = Symbol),
  #     size = 2,
  #     box.padding = unit(0.35, "lines"),
  #     point.padding = unit(0.3, "lines")
  #   )
  # ggsave(filename =file.path(resDir, "volcano_text_green.jpeg"),dpi = 300)
  
  volcano_text_blue <- volcano_2 +
    geom_text_repel(
      data = subset(dat4_dep, Accession %in% text$Accession ),
      aes(label = Symbol),
      size = 2,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  ggsave(filename =file.path(resDir, "volcano_text_blue.jpeg"),dpi = 300)
}else{stop("we can not find the file:Volcano_label_text.xlsx")}

