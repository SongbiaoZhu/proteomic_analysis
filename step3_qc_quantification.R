# step3_qc_quantification

rm(list = ls())

source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")

library(ggplot2)

# How many proteins has quantative ratio?
sum(is.na(dat1_symbol$ratio129.126 ))
sum(is.na(dat1_symbol$ratio130.127 ))
sum(is.na(dat1_symbol$ratio131.128 ))
ratio.na <- is.na(
  dat1_symbol[,c("ratio129.126","ratio130.127","ratio131.128")])
table(rowSums(ratio.na))[2]
print(paste0(table(rowSums(ratio.na))[1]," proteins have ratio" ))
print(paste0(table(rowSums(ratio.na))[2]," proteins don't have ratio" ))

print(paste0(dim(dat3_omna)[1]," proteins with ratios have more than 2 unique peptides" ))

print(paste0(dim(dat2_uni2)[1]-dim(dat3_omna)[1], " Proteins with more than 2 unique peptides but don't have ratio"))
dat2_uni2[which(is.na(dat2_uni2$ratio129.126)),][,c("Accession","Description","Unique.Peptides","ratio129.126","ratio130.127","ratio131.128","Score")]


# repeatability plot for QC
theme_set(theme_gray()+theme(axis.line = element_line(size=0.5),panel.background = element_rect(fill=NA,size=rel(20)), panel.grid.minor = element_line(colour = NA), axis.text = element_text(size=16), axis.title = element_text(size=18)))

reprodu_1 <- ggplotRegression(lm(ratio129.126 ~ ratio130.127, data = dat2_uni2[, 3:5]))+
  labs(x = "Replicate 1", y = "Replicate 2")
ggsave("out_folder/reprodu_1.jpeg",
       dpi = 600,
       width = 8.5,
       height = 8.5)
reprodu_2 <-
  ggplotRegression(lm(ratio131.128 ~ ratio130.127, data = dat2_uni2[, 3:5])) +
  labs(x = "Replicate 3", y = "Replicate 2")
ggsave("out_folder/reprodu_2.jpeg",
       dpi = 600,
       width = 8.5,
       height = 8.5)
reprodu_3 <-
  ggplotRegression(lm(ratio129.126 ~ ratio131.128, data = dat2_uni2[, 3:5])) +
  labs(x = "Replicate 1", y = "Replicate 3")
ggsave("out_folder/reprodu_3.jpeg",
       dpi = 600,
       width = 8.5,
       height = 8.5)

library(GGally)
reprodu_pairs <- ggpairs(log2(na.omit(dat2_uni2[, 3:5])),
                         upper = list(continuous = wrap("cor", size = 8)),
                         diag = list(continuous = "densityDiag"),
                         lower = list(continuous = "smooth"),
                         columnLabels = c("replicate_1", "replicate_2", "replicate_3"))
ggsave(
  "out_folder/reprodu_pairs.jpeg",
  reprodu_pairs,
  dpi = 600,
  width = 8.5,
  height = 8.5
)
reprodu_pairs_internal <- ggpairs(log2(na.omit(dat2_uni2[, 3:5])),
                         upper = list(continuous = wrap("cor", size = 8)),
                         lower = list(continuous = "smooth"),
                         axisLabels="internal")
ggsave(
  "out_folder/reprodu_pairs_internal.jpeg",
  reprodu_pairs_internal,
  dpi = 600,
  width = 8.5,
  height = 8.5
)

# PCA plot
library("FactoMineR")
library("factoextra") 
dat_pca <- t(dat3_omna[,paste0("scaledAbun",126:131)])
dat_pca <- as.data.frame(dat_pca)
group_list <- c(paste0("Control-",c(1:3)),paste0("VHL-OE-",c(1:3)))
dat_pca <- cbind(dat_pca,group_list) 
# before PCA analysis
result.pca <- PCA(dat_pca[,-ncol(dat_pca)], graph = FALSE)#现在dat_pca最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(result.pca,
             geom.ind =  "point", # show points only (nbut not "text")
             col.ind = dat_pca$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "sample"
)
ggsave("out_folder/PCA.jpeg",dpi = 300)

# heatmap plot cluster
# use the top 1000 SD proteins
cg=names(tail(sort(apply(dat3_omna[,paste0("scaledAbun",126:131)],1,sd)),1000))
library(pheatmap)
pheatmap(dat3_omna[,paste0("scaledAbun",126:131)][cg,],show_colnames =F,show_rownames = F) 
n=t(scale(t(dat3_omna[,paste0("scaledAbun",126:131)][cg,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group=c(rep("Control",3),rep("VHL-OE",3)))
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =F,show_rownames = F,cellwidth = 60,
         annotation_col=ac,filename = 'out_folder/heatmap_top1000_sd.jpeg')
