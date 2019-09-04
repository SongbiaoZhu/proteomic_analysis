# step9_gsea
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

suppressMessages(library(GSEABase))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
library(enrichplot)
library(ggplot2)
library(openxlsx)

# GSEA analysiss
# 1.read data
temp <- dat4_dep[, c("Accession", "log2FC")]
# 2. uniprot ID to Entrez_id
uniprot_id <- temp$Accession
entrez_id = bitr(uniprot_id,
                 fromType = "UNIPROT",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
head(entrez_id)
entrez_id$duplicate <- duplicated(entrez_id$ENTREZID)
summary(entrez_id$duplicate)
dat5_gsea <- merge(temp, entrez_id, by.x = "Accession", by.y = "UNIPROT")
# for simple, delete the duplicate
dat5_gsea <-  subset(dat5_gsea, duplicate == FALSE)
head(dat5_gsea, 2)
# 3. format the genelist for GSEA analysis, need ID column and FC column
## feature 1: numeric vector
geneList = dat5_gsea[, 2]
## feature 2: named vector
names(geneList) = as.character(dat5_gsea[, 3])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

# 4. read .gmt, downloaded MSigDB gene set file
if(file.exists(file.path(publicDir,"c5.all.v6.2.entrez.gmt"))){
  gmtfile <- file.path(publicDir,"c5.all.v6.2.entrez.gmt")
  c5 <- read.gmt(gmtfile)
}else{stop("we can not find the file:c5.all.v6.2.entrez.gmt")}

# 5. run GSEA
gsea_c5 <-
  GSEA(
    geneList,
    TERM2GENE = c5,
    verbose = FALSE,
    pvalueCutoff = 0.05
  )
# order GSEA result by enrichmentScore,save data
gsea_c5_sort <-
  gsea_c5[order(gsea_c5$enrichmentScore, decreasing = T)]
head(gsea_c5_sort$Description)
head(gsea_c5_sort$ID)
dim(gsea_c5_sort)
gsea_c5_sort$ID<- unlist(lapply(gsea_c5_sort$ID, function(x)firstup(x)))
gsea_c5_sort$Description<- unlist(lapply(gsea_c5_sort$Description, function(x)firstup(x)))
write.xlsx(gsea_c5_sort, file = file.path(resDir, "gsea_c5_sort.xlsx"), row.names = FALSE)
# 6. Visualization, save plot
### 将ES排名前3的画在一张图上
gseaplot2(gsea_c5, row.names(gsea_c5_sort)[1:3], pvalue_table = F)
ggsave(filename = file.path(resDir, "gsea_top3_signature.jpeg"))
# save data
save(dat5_gsea, geneList, gsea_c5,gsea_c5_sort, file = "output_step9.RData")
