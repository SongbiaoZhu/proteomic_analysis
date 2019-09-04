# step8_gopathway
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
library(clusterProfiler)
library(org.Hs.eg.db)
temp <- bitr(unique(dat4_dep$Symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
head(temp)
DEG<- dat4_dep
head(DEG)
DEG <- merge(DEG,temp,by.x='Symbol', by.y='SYMBOL')
head(DEG)
colnames(DEG)

gene_up= DEG[DEG$regulation == 1,'ENTREZID'] 
gene_down=DEG[DEG$regulation == -1,'ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )

boxplot(DEG$log2FC)

geneList=DEG$log2FC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)

## KEGG pathway analysis
if(T){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  dotplot(kk.up )
  ggsave(file.path(resDir,"kk.up.dotplot.jpeg"))
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  dotplot(kk.down )
  ggsave(file.path(resDir,"kk.down.dotplot.jpeg"))
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  dotplot(kk.diff )
  ggsave(file.path(resDir,"kk.diff.dotplot.jpeg"))
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = file.path(resDir,"kegg_up_down.jpeg"))
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = file.path(resDir,"kegg_up_down_gsea.jpeg"))
  
}        