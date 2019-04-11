# step8_gsea
rm(list=ls())
source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")
load("output_step4.RData")

suppressMessages(library(GSEABase))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
library(enrichplot)
library(ggplot2)

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
if(file.exists("c5.all.v6.2.entrez.gmt")){
  gmtfile <- "c5.all.v6.2.entrez.gmt"
  c5 <- read.gmt(gmtfile)
}else{stop("we can not find the file:c5.all.v6.2.entrez.gmt")}

# 5. run GSEA
gsea_c5 <-
  GSEA(
    geneList,
    TERM2GENE = c5,
    verbose = FALSE,
    pvalueCutoff = 0.1
  )
# order GSEA result by enrichmentScore,save data
gsea_c5_sort <-
  gsea_c5[order(gsea_c5$enrichmentScore, decreasing = T)]
head(gsea_c5_sort$Description)
head(gsea_c5_sort$ID)
dim(gsea_c5_sort)
openxlsx::write.xlsx(gsea_c5_sort, "out_folder/gsea_c5_sort.xlsx", row.names = FALSE)
# 6. Visualization, save plot
### 将ES排名前3的画在一张图上
gseaplot2(gsea_c5, row.names(gsea_c5_sort)[1:3], pvalue_table = TRUE)
ggsave("out_folder/gsea_top3_signature.pdf")
### 将某一类，eg: IFN related 画在一张图上
(index_1 <- which(grepl("INTERFERON", gsea_c5_sort$ID)))
gsea_c5_sort$ID[index_1]
gseaplot2(
  gsea_c5,
  row.names(gsea_c5_sort)[index_1[c(1, 5)]],
  pvalue_table = FALSE,
  color = c("red", "green"),
  base_size = 14,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  ES_geom = "line"
)
ggsave("out_folder/gsea_ifn.pdf")
### 将APC related 画在一张图上
(index_2 <- which(grepl("ANTIGEN_PROCESSING", gsea_c5_sort$ID)))
gsea_c5_sort$ID[index_2]
gseaplot2(
  gsea_c5,
  row.names(gsea_c5_sort)[index_2[c(1, 2)]],
  pvalue_table = FALSE,
  color = c("red", "green"),
  base_size = 14,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  ES_geom = "line"
)
ggsave("out_folder/gsea_apc.pdf")

# save data
save(dat5_gsea, geneList, gmtfile, gsea_c5,gsea_c5_sort, file = "output_step8.RData")
