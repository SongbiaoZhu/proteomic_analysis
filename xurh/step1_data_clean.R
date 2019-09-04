rm(list = ls())
options(stringsAsFactors = F)

source("udf.R")
creDir("data")
creDir("res")
creDir("public")

dataDir = file.path(getwd(),'data')
resDir = file.path(getwd(),'res')
publicDir = file.path(getwd(),'public')

library(stringr)
library(openxlsx)
library(dplyr)
# read data txt
dat <-
  read.xlsx(file.path(dataDir, 'MS171081_17008、17012、17013.xlsx'),sheet = 1, 
            colNames = TRUE,na.strings = "NA")

# count how many uncharacterized proteins, but not remove
dat_unchar <- dat %>%
  mutate(unchar = grepl("Uncharacterized protein", Description)) %>%
  filter(unchar == TRUE)
head(dat_unchar$Description,10)

# extract key columns, add gene symbol
dat1_symbol <- dat[, c(1, 2, 6, 9, 12, 15, 3, 4)]#提取关键列
colnames(dat1_symbol) <- c(
  "Accession",
  "Description",
  "Unique.Peptides" ,
  paste0("ratio", c("127.126", "129.128", "131.130")),
  "Score",
  "Coverage"
)
dat1_symbol$log2ratio127.126 = log2(dat1_symbol$ratio127.126)
dat1_symbol$log2ratio129.128 = log2(dat1_symbol$ratio129.128)
dat1_symbol$log2ratio131.130 = log2(dat1_symbol$ratio131.130)

temp <- str_split_fixed(dat1_symbol$Description, "GN=", 2)
temp <- str_split_fixed(temp[, 2], " ", 2)
temp <- str_split_fixed(temp[, 1], "-", 2)
grep("_HUMAN", dat1_symbol$Symbol)
dat1_symbol <- dat1_symbol %>%
  mutate(Symbol = temp[, 1])
sum(duplicated(dat1_symbol$Symbol))

# filter unipep>=2
dat2_uni2 <- dat1_symbol %>%
  filter(Unique.Peptides >= 2)

# remove na in the ratio columns
apply(dat2_uni2, 2, function(x)sum(is.na(x)))
dat3_omna <- na.omit(dat2_uni2)

# save cleaned data
save(dat,dat_unchar,dat1_symbol, dat2_uni2, dat3_omna, file = "output_step1.RData")
