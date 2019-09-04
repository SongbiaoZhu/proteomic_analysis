rm(list=ls())

# edit the step1.R, change the dataset name
scripts <- list.files(pattern = "^step")
sapply(scripts, source) 
