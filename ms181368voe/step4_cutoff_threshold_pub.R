
rm(list=ls())

source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")

# calculate ratio variance to see whether change the fc cut
library(matrixStats)
library(dplyr)

dat3_omna_var <- dat3_omna 
dat3_omna_var $abun.SD <- rowSds(as.matrix(dat3_omna[, c(9:14)]), na.rm = TRUE)
dat3_omna_var $abun.mean <- 1/6 * (
  dat3_omna_var[,9] + dat3_omna_var[,10] + dat3_omna_var[,11] +
    dat3_omna_var[,12] + dat3_omna_var[,13] + dat3_omna_var[,14] )
dat3_omna_var$rel.var <- dat3_omna_var$abun.SD / dat3_omna_var$abun.mean

rel.var.table <- as.data.frame(table(cut(
  dat3_omna_var$rel.var ,
  breaks = seq(0, 1, 0.1),
  include.lowest = TRUE
)))
rel.var.table$CumFreq <- cumsum(rel.var.table$Freq)
rel.var.table$percent <- rel.var.table$Freq/sum( rel.var.table$Freq)
rel.var.table$Cumulative <- cumsum(rel.var.table$percent)

library(ggplot2)

gg <-
  ggplot(data = rel.var.table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.8) +
  xlab("Relative variation") +
  ylab("Protein counts") 
#plot line on same graph
# rate multiplied by 10000000 to get on same scale as bars
gg <-
  gg +  geom_line(data = rel.var.table,
                  aes(x = Var1,
                      y = (Cumulative) * sum(rel.var.table$Freq),
                             group = 1),
                  color="black",size= 0.8, inherit.aes = FALSE) +
  scale_y_continuous(expand = c(0, 0),
                      sec.axis = sec_axis( ~ . /  sum(rel.var.table$Freq), breaks=c(0.25,0.5,0.75,1.0, 1.25),name = "Cumulative coverage")) +
  geom_point(data = rel.var.table,
             aes(
               x = Var1,
               y = (Cumulative) *  sum(rel.var.table$Freq),
               group = 1
             ),
             inherit.aes = FALSE)
gg
library(ggThemeAssist)
# ggThemeAssistGadget(gg)

gg + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(linetype = "solid", size = 1), 
    axis.ticks = element_line(colour = "gray0"), 
    axis.title = element_text(size = 12, 
        face = "bold"), axis.text = element_text(size = 12), 
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, face = "bold"), 
    axis.text.y = element_text(size = 12, face = "bold"), 
    panel.background = element_rect(fill = NA), 
    legend.position = "none") +labs(size = 12)


ggsave("out_folder/rel.var_bar_pub.jpeg")
