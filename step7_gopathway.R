# step7_gopathway
rm(list=ls())
source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")
load("output_step4.RData")

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
  dotplot(kk.up );ggsave('out_folder/kk.up.dotplot.png')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  dotplot(kk.down );ggsave('out_folder/kk.down.dotplot.png')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  dotplot(kk.diff );ggsave('out_folder/kk.diff.dotplot.png')
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  source('functions.R')
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = 'out_folder/kegg_up_down.png')
  
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
  ggsave(g_kegg,filename = 'out_folder/kegg_up_down_gsea.png')
  
  
}

### GO database analysis 
{
  
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  
  if(F){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file = 'go_enrich_results.Rdata')
    
  }
  
  
  load(file = 'go_enrich_results.Rdata')
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0('out_folder/dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
  
}
save(DEG,kk.diff,kk.down,kk.up,kk_gse,file = 'output_step7.Rdata')
