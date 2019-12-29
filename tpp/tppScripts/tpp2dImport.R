tpp2dImport
function (configTable = NULL, data = NULL, idVar = "gene_name", 
    addCol = NULL, intensityStr = "signal_sum_", qualColName = "qupm", 
    nonZeroCols = "qssm", fcStr = NULL) 
{
    configWide <- importCheckConfigTable(infoTable = configTable, 
        type = "2D")
    experiment = unique_ID <- NULL
    message("Importing data...")
    dataList <- import2DTR_main(configTable = configWide, data = data, 
        idVar = idVar, fcStr = fcStr, addCol = addCol, naStrs = c("NA", 
            "n/d", "NaN"), intensityStr = intensityStr, qualColName = qualColName, 
        nonZeroCols = nonZeroCols)
    dataWide <- importFct_createCCRInputFrom2DData(configTable = configWide, 
        data.list = dataList, intensityStr = intensityStr, fcStr = fcStr) %>% 
        mutate(experiment = factor(experiment), unique_ID = factor(unique_ID))
    dataWide <- dataWide %>% arrange_(.dots = c(idVar, "temperature"))
    attr(dataWide, "configTable") <- configWide
    importSettings <- list(proteinIdCol = idVar, uniqueIdCol = "unique_ID", 
        addCol = addCol, intensityStr = intensityStr, qualColName = qualColName, 
        nonZeroCols = nonZeroCols, fcStr = fcStr)
    attr(dataWide, "importSettings") <- importSettings
    message("Done.")
    return(dataWide)
}
<bytecode: 0x00000000101ad688>
<environment: namespace:TPP>