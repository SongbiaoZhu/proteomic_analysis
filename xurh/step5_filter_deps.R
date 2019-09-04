
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
# filter DEPs
fc_cut <- 1.5
p_cut <- 0.05
dat4_dep <- dat3_omna
dat4_dep$mean.ratio <- apply(dat4_dep[,c(4:6)],1,mean)
dat4_dep$log2FC <- log2(dat4_dep$mean.ratio)
dat4_dep$p.value <-
  apply(dat4_dep[,c(9:11)], 1, function(x)
    t.test(x[1:3], mu = 0)$p.value)
dat4_dep$ad.p.value <- p.adjust(dat4_dep$p.value, "fdr")
#The "BH" (aka "fdr") , control the false discovery rate
if(F){
  dat4_dep$threshold = as.factor(dat4_dep$ad.p.value < p_cut &
                                   (dat4_dep$log2FC > log2(fc_cut) |
                                      (dat4_dep$log2FC < log2(1 / fc_cut))))
  dat4_dep <- within(dat4_dep, {
    regulation = as.factor (ifelse(threshold == FALSE, 0,
                                   ifelse(
                                     log2FC < log2(1 / fc_cut), -1,
                                     ifelse(log2FC > log2(fc_cut), 1, 0)
                                   )))})
}
dat4_dep$threshold = as.factor(dat4_dep$p.value < p_cut &
                                 (dat4_dep$log2FC > log2(fc_cut) |
                                    (dat4_dep$log2FC < log2(1 / fc_cut))))
dat4_dep <- within(dat4_dep, {
  regulation = as.factor (ifelse(threshold == FALSE, 0,
                                 ifelse(
                                   log2FC < log2(1 / fc_cut), -1,
                                   ifelse(log2FC > log2(fc_cut), 1, 0)
                                 )))})
# save filtered data
save(dat4_dep,fc_cut,p_cut, file = "output_step5.RData")

