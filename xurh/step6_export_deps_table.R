# step6_export_deps_tables
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

library(openxlsx)
library(dplyr)
library(stringr)

S1.Table <- dat4_dep %>%
  filter(as.character(regulation) ==  1) %>%
  mutate("Unique_Peptides" = Unique.Peptides) %>%
  mutate("Description" = str_split_fixed(Description, " OS=", 2)[,1]) %>%
  mutate("Average_Ratio" = mean.ratio) %>%
  mutate("p.value" = p.value) %>%
  select ("Accession",
          "Description",
          "Score",
          "Coverage",
          "Unique_Peptides",
          "Average_Ratio",
          "p.value")  %>%
  arrange(desc(Average_Ratio))

S2.Table <- dat4_dep %>%
  filter(as.character(regulation) ==  -1) %>%
  mutate("Unique_Peptides" = Unique.Peptides) %>%
  mutate("Description" = str_split_fixed(Description, " OS=", 2)[,1]) %>%
  mutate("Average_Ratio" = mean.ratio) %>%
  mutate("p.value" = p.value) %>%
  select ("Accession",
          "Description",
          "Score",
          "Coverage",
          "Unique_Peptides",
          "Average_Ratio",
          "p.value")  %>%
  arrange(Average_Ratio)

deps <- dat4_dep %>%
  filter(regulation ==
           1 | regulation == -1) %>%
  arrange(desc(mean.ratio))
deps_up <- dat4_dep %>%
  filter(regulation ==  1) %>%
  arrange(desc(mean.ratio))
deps_down <- dat4_dep %>%
  filter(regulation ==  -1) %>%
  arrange(mean.ratio)

write.xlsx(dat, file.path(resDir,"dat_Proteins.xlsx"), row.names = FALSE)
write.xlsx(dat1_symbol, file.path(resDir,"dat1_symbol.xlsx"), row.names = FALSE)
write.xlsx(deps,file.path(resDir,"deps_proteins.xlsx"),row.names = FALSE)
write.xlsx(deps_up,file.path(resDir,"deps_up_proteins.xlsx"),row.names = FALSE)
write.xlsx(deps_down,file.path(resDir,"deps_down_proteins.xlsx"),row.names = FALSE)
write.xlsx(dat4_dep[, c("Accession", "log2FC", "p.value")], file.path(resDir,"forIPA.xlsx"), row.names = FALSE)
write.xlsx(S1.Table, file.path(resDir,"S1.Table.Up-regulated.xlsx"), row.names = FALSE)
write.xlsx(S2.Table, file.path(resDir,"S2.Table.Down-regulated.xlsx"), row.names = FALSE)

write.table(
  deps_up[, "Accession"],
  file.path(resDir,"deps_up_accession.txt"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  deps_down[, "Accession"],
  file.path(resDir,"deps_down_accession.txt"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  deps_up [, "Symbol"],
  file.path(resDir,"deps_up_Symbol.txt"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  deps_down [, "Symbol"],
  file.path(resDir,"deps_down_Symbol.txt"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  subset(dat4_dep[, c("Symbol", "ratio127.126", "ratio129.128", "ratio131.130")],!is.na(Symbol)),
  file.path(resDir,"forPSEA.txt"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
