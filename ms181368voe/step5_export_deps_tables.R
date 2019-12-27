# step5_export_deps_tables
rm(list=ls())

source("functions.R")
cre_folder("out_folder")

load("output_step1.RData")
load("output_step4.RData")

library(openxlsx)
library(dplyr)
library(stringr)

S1.Table <- dat4_dep %>%
  filter(as.character(regulation) ==  1) %>%
  mutate("Unique_Peptides" = Unique.Peptides) %>%
  mutate("Description" = str_split_fixed(Description, " OS=", 2)[,1]) %>%
  mutate("Average_Ratio" = mean.ratio) %>%
  mutate("adjusted_p_value" = ad.p.value) %>%
  select ("Accession",
          "Description",
          "Score",
          "Coverage",
          "Unique_Peptides",
          "Average_Ratio",
          "adjusted_p_value")  %>%
  arrange(desc(Average_Ratio))

S2.Table <- dat4_dep %>%
  filter(as.character(regulation) ==  -1) %>%
  mutate("Unique_Peptides" = Unique.Peptides) %>%
  mutate("Description" = str_split_fixed(Description, " OS=", 2)[,1]) %>%
  mutate("Average_Ratio" = mean.ratio) %>%
  mutate("adjusted_p_value" = ad.p.value) %>%
  select ("Accession",
          "Description",
          "Score",
          "Coverage",
          "Unique_Peptides",
          "Average_Ratio",
          "adjusted_p_value")  %>%
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

write.xlsx(dat, "out_folder/dat_Proteins_81cols.xlsx", row.names = FALSE)
write.xlsx(dat1_symbol, "out_folder/dat1_symbol.xlsx", row.names = FALSE)
write.xlsx(dat4_dep, "out_folder/dat4_dep.xlsx", row.names = FALSE)

write.xlsx(deps,
           "out_folder/deps_proteins.xlsx",
           row.names = FALSE)
write.xlsx(deps_up,
           "out_folder/deps_up_proteins.xlsx",
           row.names = FALSE)
write.xlsx(deps_down,
           "out_folder/deps_down_proteins.xlsx",
           row.names = FALSE)

write.xlsx(dat4_dep[, c("Accession", "log2FC", "ad.p.value")], "out_folder/forIPA.xlsx", row.names = FALSE)
write.xlsx(S1.Table, "out_folder/S1.Table.Up-regulated.xlsx", row.names = FALSE)
write.xlsx(S2.Table, "out_folder/S2.Table.Down-regulated.xlsx", row.names = FALSE)

write.table(
  deps_up[, "Accession"],
  "out_folder/deps_up_accession.txt",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  deps_down[, "Accession"],
  "out_folder/deps_down_accession.txt",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  deps_up [, "Symbol"],
  "out_folder/deps_up_Symbol.txt",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  deps_down [, "Symbol"],
  "out_folder/deps_down_Symbol.txt",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
write.table(
  subset(dat4_dep[, c("Symbol", "ratio129.126", "ratio130.127", "ratio131.128")],!is.na(Symbol)),
  "out_folder/forPSEA.txt",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
