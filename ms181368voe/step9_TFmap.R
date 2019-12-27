# setp8_TFmap
rm(list=ls())
source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")
load("output_step4.RData")

# 1st, import the TF database
if(file.exists("TF_DatabaseExtract_v_1.01.csv")){
hTFs <- read.csv("TF_DatabaseExtract_v_1.01.csv", sep = ",", stringsAsFactors = FALSE)
hTFs <- subset(hTFs, Is.TF. == "Yes")
}else{stop("we can not find the file:TF_DatabaseExtract_v_1.01.csv")}
# 2, see how many identified TFs, how many upTFs AND downTFs
tf.dat <- merge(hTFs, dat1_symbol, by.x = "HGNC.symbol", by.y = "Symbol")
tf.deps <-
  merge(hTFs, dat4_dep[dat4_dep$threshold == TRUE, ], by.x = "HGNC.symbol", by.y = "Symbol")
tf.updeps <-
  merge(hTFs, dat4_dep[dat4_dep$regulation == 1, ], by.x = "HGNC.symbol", by.y = "Symbol")
tf.downdeps <-
  merge(hTFs, dat4_dep[dat4_dep$regulation == -1, ], by.x = "HGNC.symbol", by.y = "Symbol")
# 3, export tf tables
library(openxlsx)
write.xlsx(tf.dat, "out_folder/Identified_TFs.xlsx", row.names = FALSE)
write.xlsx(tf.deps, "out_folder/deps_TFs.xlsx", row.names = FALSE)
write.xlsx(tf.updeps, "out_folder/updeps_TFs.xlsx", row.names = FALSE)
write.xlsx(tf.downdeps, "out_folder/downdeps_TFs.xlsx", row.names = FALSE)
# 4,venn plot
suppressMessages(if(!require(VennDiagram)){
  install.packages("VennDiagram")
  library(VennDiagram)
})
## how many identified tfs
venn.diagram(
  x = list(dat1_symbol$Symbol, hTFs$HGNC.symbol ),
  category.names = c("Proteomics" , "hTFs " ),
  filename = 'out_folder/Venn_TFs_Identified.tiff',
  na = "remove",
  output = TRUE ,
  imagetype = "tiff",
  height = 1000,
  width = 1000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = c('solid', 'longdash'),
  cex = 0.8,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.0,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.default.pos = "text",
  cat.pos = c(-5, 5)
)
## how many identified deps tfs
venn.diagram(
  x = list(dat4_dep[dat4_dep$threshold == TRUE,]$Symbol, hTFs$HGNC.symbol ),
  category.names = c("deps" , "hTFs " ),
  filename = 'out_folder/Venn_TFs_deps.tiff',
  na = "remove",
  output = TRUE ,
  imagetype = "tiff",
  height = 1000,
  width = 1000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = c('solid', 'longdash'),
  cex = 0.8,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.0,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.default.pos = "text",
  cat.pos = c(-5, 5)
)
# Quadruple Venn plots 
quad.venn <- list(
  All = dat1_symbol$Symbol,
  hTFs = hTFs$HGNC.symbol,
  upTFs = tf.updeps$HGNC.symbol,
  downTFs = tf.downdeps$HGNC.symbol
)
venn.diagram(
  quad.venn,
  filename = "out_folder/Venn_Quad.tiff",
  height = 3000,
  width = 3000,
  na = "remove",
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)

# heatmap plot cluster
# use the tf.deps
library(pheatmap)
tf.deps_heat <- tf.deps[,c("HGNC.symbol",paste0("scaledAbun",126:131))]
rownames(tf.deps_heat) <- tf.deps_heat[,1]
tf.deps_heat <- tf.deps_heat[,-1]
pheatmap(tf.deps_heat,show_colnames =F,show_rownames = F) 
n=t(scale(t(tf.deps_heat))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group=c(rep("Control",3),rep("VHL-OE",3)))
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =F,show_rownames = T,display_numbers = F,fontsize_row = 10, annotation_legend = T,
         cellwidth = 70,cellheight = 20,
         annotation_col=ac,filename = 'out_folder/heatmap_tf_deps.jpeg')

# save data
save(hTFs,tf.dat,tf.deps, tf.updeps, tf.downdeps, tf.deps_heat, file = "output_step9.RData")


