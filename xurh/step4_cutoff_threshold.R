# step4_cutoff_threshold
rm(list = ls())
options(stringsAsFactors = F)
source("udf.R")
creDir("data")
creDir("res")
creDir("public")
dataDir = file.path(getwd(),'data')
resDir = file.path(getwd(),'res')
publicDir = file.path(getwd(),'public')
library(ggplot2)
library(ggThemeAssist)
load("output_step1.RData")

# calculate ratio variance to see whether change the fc cut
library(matrixStats)
library(dplyr)
dat3_omna_var <- dat3_omna 
dat3_omna_var $ratio.SD <- rowSds(as.matrix(dat3_omna[, c(4:6)]), na.rm = TRUE)
dat3_omna_var $ratio.mean <- 1/3 * (
  dat3_omna_var[,4] + dat3_omna_var[,5] + dat3_omna_var[,6])
dat3_omna_var$rel.var <- dat3_omna_var$ratio.SD / dat3_omna_var$ratio.mean

rel.var.table <- as.data.frame(table(cut(
  dat3_omna_var$rel.var ,
  breaks = seq(0, 1, 0.1),
  include.lowest = TRUE
)))
rel.var.table$CumFreq <- cumsum(rel.var.table$Freq)
rel.var.table$percent <- rel.var.table$Freq/sum( rel.var.table$Freq)
rel.var.table$Cumulative <- cumsum(rel.var.table$percent)

# plot
library(ggplot2)
if(T){
  gg <-
    ggplot(data = rel.var.table, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", width = 0.8) + 
    xlab("Relative variation") +
    ylab("Protein counts") 
  gg
  
  gg <- gg + geom_point(data = rel.var.table,
                        aes(
                          x = Var1,
                          y = (Cumulative) *  sum(rel.var.table$Freq),
                          group = 1
                        ), size = 0.9, color="black", 
                        inherit.aes = TRUE) 
  gg
  
  gg <- gg + geom_line(data = rel.var.table,
                       aes(x = Var1,
                           y = (Cumulative) * sum(rel.var.table$Freq),
                           group = 1), size = 0.7, color="black", 
                       inherit.aes = TRUE) 
  gg
  
  gg <- gg + scale_y_continuous(
    sec.axis = sec_axis( ~ . /  sum(rel.var.table$Freq),
                         name = "Cumulative coverage"))
  gg
  
  gg <- gg + theme(plot.subtitle = element_text(vjust = 1), 
                   plot.caption = element_text(vjust = 1), 
                   axis.line = element_line(linetype = "solid", size = 1), 
                   axis.ticks = element_line(colour = "gray0"), 
                   axis.title = element_text(size = 12, 
                                             face = "bold"), axis.text = element_text(size = 12), 
                   axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, face = "bold"), 
                   axis.text.y = element_text(size = 12, face = "bold"), 
                   panel.background = element_rect(fill = NA), 
                   legend.position = "none") +labs(size = 12)
  gg
  ggsave(filename = file.path(resDir, "rel.var_bar_pub.jpeg") , dpi = 300)
  
}
