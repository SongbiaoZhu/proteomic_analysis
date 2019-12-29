tpp2dNormalize
function (configTable = NULL, data, fcStr = NULL) 
{
    if (!missing(configTable)) {
        warning("`configTable` is deprecated.", call. = TRUE)
    }
    if (!missing(fcStr)) {
        warning("`fcStr` is deprecated.", call. = TRUE)
    }
    median_per_temp_and_conc = y = yNormalized = columnName <- NULL
    checkFunctionArgs(match.call(), c("data"))
    importSettings <- attr(data, "importSettings")
    idVar <- checkAndReturnDataSetting(importSettings, "proteinIdCol", 
        colnames(data))
    fcStr <- checkAndReturnDataSetting(importSettings, "fcStr", 
        colnames(data))
    fcStrNorm <- paste("norm", fcStr, sep = "_")
    fcCols <- grep(fcStr, colnames(data), value = TRUE)
    fcNormCols <- gsub(fcStr, fcStrNorm, fcCols)
    message("Performing median normalization per temperature...")
    dataLong <- data %>% subset(select = c(idVar, "temperature", 
        fcCols)) %>% gather_("columnName", "y", gather_cols = fcCols) %>% 
        mutate(y = as.numeric(as.character(y)))
    normCoeffs <- dataLong %>% group_by_(.dots = c("temperature", 
        "columnName")) %>% summarise(median_per_temp_and_conc = median(y, 
        na.rm = TRUE))
    dataNormed <- dataLong %>% left_join(normCoeffs, by = c("temperature", 
        "columnName")) %>% mutate(yNormalized = y/median_per_temp_and_conc) %>% 
        select(-median_per_temp_and_conc, -y) %>% mutate(columnName = gsub(fcStr, 
        fcStrNorm, columnName) %>% factor(levels = fcNormCols)) %>% 
        spread(columnName, yNormalized)
    out <- data %>% left_join(dataNormed, by = c(idVar, "temperature"))
    rownames(out) <- attr(data, "row.names")
    importSettings$fcStrNorm = fcStrNorm
    attr(out, "importSettings") <- importSettings
    attr(out, "configTable") <- attr(data, "configTable")
    message("Done.")
    return(out)
}
<bytecode: 0x00000000173ad280>
<environment: namespace:TPP>