#read PD export data

rm(list=ls())
library(stringr)
library(dplyr)
# read data txt
dat <-
  read.table(
    "MS181368-PVOE-3-channelAver_Proteins.txt",
    header = TRUE,
    stringsAsFactors = FALSE
  )

# count how many uncharacterized proteins, but not remove
dat_unchar <- dat %>%
  mutate(unchar = grepl("Uncharacterized protein", Description)) %>%
  filter(unchar == TRUE)
head(dat_unchar$Description,10)

# extract key columns, add gene symbol
dat1_symbol <- dat[, c(6, 7, 17, 46:75, 78, 79)]#提取关键列
colnames(dat1_symbol) <- c(
  "Accession",
  "Description",
  "Unique.Peptides" ,
  paste0("ratio", c("129.126", "130.127", "131.128")),
  paste0("log2ratio", c("129.126", "130.127", "131.128")),
  paste0("scaledAbun", c(126:131)),
  paste0("normAbun", c(126:131)),
  paste0("Abun", c(126:131)),
  paste0("Abuncount", c(126:131)),
  "Score",
  "Coverage"
)
temp <- str_split_fixed(dat1_symbol$Description, "GN=", 2)
temp <- str_split_fixed(temp[, 2], "PE=", 2)
temp <- gsub(" ", "", temp)
dat1_symbol <- dat1_symbol %>%
  mutate(Symbol = temp[, 1])
sum(duplicated(dat1_symbol$Symbol))

# filter unipep>=2
dat2_uni2 <- dat1_symbol %>%
  filter(Unique.Peptides >= 2) %>%
  select (c(1,3:15, 2, 34:36))

# remove na in the ratio columns
dat2_uni2[which(is.na(dat2_uni2$ratio129.126)),][,c(3:5,15,16,17)]
dat3_omna <- na.omit(dat2_uni2)

# save cleaned data
save(dat,dat_unchar,dat1_symbol, dat2_uni2, dat3_omna, file = "output_step1.RData")


