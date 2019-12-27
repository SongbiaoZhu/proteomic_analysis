# step11_heatmap_custom
rm(list = ls())
source("functions.R")
cre_folder("out_folder")

library(ggplot2)
theme_set(
  theme_gray() + theme(
    axis.line = element_line(size = 0.5),
    panel.background = element_rect(fill = NA, size = rel(20)),
    panel.grid.minor = element_line(colour = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )
)

if (file.exists("heatmap_custom.xlsx")) {
  dat_process <- openxlsx::read.xlsx("heatmap_custom.xlsx", sheet = 1)
  library(pheatmap)
  dat_process_heat <-
    dat_process[, c("Symbol", paste0("scaledAbun", 126:131))]
  rownames(dat_process_heat) <- dat_process_heat[, 1]
  dat_process_heat <- dat_process_heat[, -1]
  pheatmap(dat_process_heat,
           show_colnames = F,
           show_rownames = F)
  n = t(scale(t(dat_process_heat)))
  n[n > 2] = 2
  n[n < -2] = -2
  n[1:4, 1:4]
  pheatmap(n, show_colnames = F, show_rownames = F)
  ac = data.frame(group = c(rep("Control", 3), rep("VHL-OE", 3)))
  rownames(ac) = colnames(n)
  pheatmap(
    n,
    show_colnames = F,
    show_rownames = T,
    display_numbers = F,
    fontsize_row = 10,
    annotation_legend = T,
    cellwidth = 70,
    cellheight = 20,
    annotation_col = ac,
    filename = 'out_folder/heatmap_process.jpeg'
  )
} else
{
  stop("we can not find the file:heatmap_custom.xlsx")
}
