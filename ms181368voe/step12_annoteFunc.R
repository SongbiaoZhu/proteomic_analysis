rm(list = ls())

source("functions.R")
cre_folder("out_folder")
library(openxlsx)

load("output_step1.RData")
load('E:/publicData/AccSymProFuncHs.RData')
datFunc <- merge(dat, AccSymProFunc, by = 'Accession')
symFreq <-table(datFunc$'Symbol')
symFreq['Nosymbol']
funcFreq <-table(datFunc$'Function')
funcFreq['Not found function summary.']

load("output_step4.RData")
deps <- dat4_dep %>%
  select(-Symbol) %>%
  filter(threshold == TRUE) %>%
  left_join(.,AccSymProFunc, by = 'Accession') %>%
  arrange(desc(mean.ratio))

deps_up <- deps %>%
  filter(regulation == 1) %>%
  arrange(desc(mean.ratio))

deps_down <- deps %>%
  filter(regulation == -1) %>%
  arrange(mean.ratio)

write.xlsx(deps,
           "out_folder/deps_func.xlsx",
           row.names = FALSE)
write.xlsx(deps_up,
           "out_folder/deps_up_func.xlsx",
           row.names = FALSE)
write.xlsx(deps_down,
           "out_folder/deps_down_func.xlsx",
           row.names = FALSE)