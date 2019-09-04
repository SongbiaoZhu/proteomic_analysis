# step3_qc_quantification
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

# How many proteins has quantative ratio?
sum(is.na(dat1_symbol$ratio127.126 ))
sum(is.na(dat1_symbol$ratio129.128 ))
sum(is.na(dat1_symbol$ratio131.130 ))
ratio.na <- is.na(
  dat1_symbol[,c("ratio127.126","ratio129.128","ratio131.130")])
table(rowSums(ratio.na))[2]
print(paste0(table(rowSums(ratio.na))[1]," proteins have ratio" ))
print(paste0(table(rowSums(ratio.na))[2]," proteins don't have ratio" ))
print(paste0(dim(dat3_omna)[1]," proteins with ratios have more than 2 unique peptides" ))
print(paste0(dim(dat2_uni2)[1]-dim(dat3_omna)[1], " Proteins with more than 2 unique peptides but don't have ratio"))
# repeatability plot for QC
# dat2 is a subset dataframe for repeatability
dat2 <- dat3_omna[,9:11]
colnames(dat2) <- c("x", "y", "z")
summary(lm(x ~ y, data = dat2))

axis.limit = c(floor(min(c(dat2$x,dat2$y))),ceiling(max(c(dat2$x,dat2$y))))
psingle <-
  ggplot(data = dat2, aes(x = x, y = y)) +
  geom_point() +
  scale_x_continuous(name = "Replicate 1", breaks = seq(axis.limit[1], axis.limit[2],1), labels = seq(axis.limit[1], axis.limit[2],1), limits = axis.limit) +
  scale_y_continuous(name = "Replicate 2", breaks = seq(axis.limit[1], axis.limit[2],1), labels = seq(axis.limit[1], axis.limit[2],1), limits = axis.limit) + 
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "red",
    se = FALSE
  )
# ggThemeAssistGadget(psingle)
psingle <- psingle + theme(plot.subtitle = element_text(vjust = 1), 
                           plot.caption = element_text(vjust = 1), 
                           axis.line = element_line(linetype = "solid"), 
                           axis.ticks = element_line(colour = "black", size = 0.5), 
                           axis.title = element_text(face = "plain"), 
                           axis.text = element_text(face = "bold"), 
                           plot.title = element_text(face = "bold"), 
                           panel.background = element_rect(fill = NA), 
                           legend.position = "none")
psingle
ggsave(filename = file.path(resDir, "repeat_scatter.jpeg") , dpi = 300)

library(GGally)
reprodu_pairs <- ggpairs(dat2,
                         upper = list(continuous = wrap("cor", size = 8)),
                         diag = list(continuous = "densityDiag"),
                         lower = list(continuous = "smooth"),
                         columnLabels = c("replicate_1", "replicate_2", "replicate_3"))
reprodu_pairs
ggsave(filename = file.path(resDir, "reprodu_pairs.jpeg") , dpi = 300)
