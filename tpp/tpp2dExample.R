rm(list = ls())
library(TPP)
data(panobinostat_2DTPP_smallExample)
config_tpp2d <- panobinostat_2DTPP_config
data_tpp2d <- panobinostat_2DTPP_data
tpp2dResults <- analyze2DTPP(configTable = config_tpp2d, 
                             data = data_tpp2d,
                             methods=c("doseResponse"),
                             createReport="none",
                             nCores='max',
                             idVar = "representative",
                             addCol = "clustername",
                             intensityStr = "sumionarea_protein_",
                             nonZeroCols = "qusm")
