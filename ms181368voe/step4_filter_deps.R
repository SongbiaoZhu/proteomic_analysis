# step4_filter deps
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
openxlsx::write.xlsx(rel.var.table, "out_folder/Variation_cutoff.xlsx", row.names = FALSE)

library(ggplot2)

rel.var_bar <-
  ggplot(data = rel.var.table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Relative Variation") +
  ylab("Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept =  sum(rel.var.table$Freq)*0.88, color = "red", linetype = "dashed")
#plot line on same graph
# rate multiplied by 10000000 to get on same scale as bars
rel.var_bar <-
  rel.var_bar +  geom_line(data = rel.var.table,
                           aes(
                             x = Var1,
                             y = (Cumulative) * sum(rel.var.table$Freq),
                             group = 1
                           ),color="black",size=1,
                           inherit.aes = FALSE) +
  scale_y_continuous(sec.axis = sec_axis( ~ . /  sum(rel.var.table$Freq), breaks=c(0.25,0.5,0.75,0.88,1.0),name = "Cumulative Coverage")) +
  geom_point(data = rel.var.table,
             aes(
               x = Var1,
               y = (Cumulative) *  sum(rel.var.table$Freq),
               group = 1
             ),
             inherit.aes = FALSE)
rel.var_bar
ggsave("out_folder/rel.var_bar.jpeg", width = 12, height = 10)

# filter DEPs
fc_cut <- 1.5
p_cut <- 0.05
dat4_dep <- dat3_omna[, c(1, 9:14, 3:8, 15:18,2)]
dat4_dep$meanCon <- apply(dat4_dep[, c(2:4)], 1, mean)
dat4_dep$meanTreat <- apply(dat4_dep[, c(5:7)], 1, mean)
dat4_dep$p.value <-
  apply(dat4_dep[,-c(1, 8:20)], 1, function(x)
    t.test(x[1:3], x[4:6])$p.value)
dat4_dep$mean.ratio <- dat4_dep$meanTreat / dat4_dep$meanCon
dat4_dep$log2FC <- log2(dat4_dep$mean.ratio)
dat4_dep$ad.p.value <- p.adjust(dat4_dep$p.value, "fdr")
#The "BH" (aka "fdr") , control the false discovery rate
dat4_dep$threshold = as.factor(dat4_dep$ad.p.value < p_cut &
                             (dat4_dep$log2FC > log2(fc_cut) |
                                (dat4_dep$log2FC < log2(1 / fc_cut))))
dat4_dep <- within(dat4_dep, {
  regulation = as.factor (ifelse(threshold == FALSE, 0,
                                 ifelse(
                                   log2FC < log2(1 / fc_cut), -1,
                                   ifelse(log2FC > log2(fc_cut), 1, 0)
                                 )))
})

# save filtered data
save(dat4_dep,fc_cut,p_cut, rel.var.table, file = "output_step4.RData")


