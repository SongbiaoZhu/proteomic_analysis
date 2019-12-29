tpp2dExport
function (configTable = NULL, tab, resultPath = NULL, idVar = NULL, 
    fcStr = NULL, intensityStr = NULL, outPath, addCol = NULL, 
    normalizedData = NULL, trRef = NULL, addPlotColumns = TRUE) 
{
    if (!missing(configTable)) {
        warning("`configTable` is deprecated.", call. = TRUE)
    }
    if (!missing(resultPath)) {
        warning("`resultPath` is deprecated.", call. = TRUE)
    }
    if (!missing(idVar)) {
        warning("`idVar` is deprecated.", call. = TRUE)
    }
    if (!missing(fcStr)) {
        warning("`fcStr` is deprecated.", call. = TRUE)
    }
    if (!missing(intensityStr)) {
        warning("`intensityStr` is deprecated.", call. = TRUE)
    }
    if (!missing(normalizedData)) {
        warning("`normalizedData` is deprecated.", call. = TRUE)
    }
    checkFunctionArgs(match.call(), c("tab", "outPath"))
    configTable <- attr(tab, "configTable")
    importSettings <- attr(tab, "importSettings")
    idVar <- checkAndReturnDataSetting(importSettings, "proteinIdCol", 
        colnames(tab))
    fcStr <- checkAndReturnDataSetting(importSettings, "fcStr", 
        colnames(tab))
    intensityStr <- checkAndReturnDataSetting(importSettings, 
        "intensityStr", colnames(tab))
    normalizedData <- !is.null(importSettings$fcStrNorm)
    fTmp <- paste0(format(Sys.time(), "%Y-%m-%d"), "_results_2D_TPP.xlsx")
    fileName <- file.path(outPath, fTmp)
    message("Writing results to file: ", fileName)
    tab <- exportFct_convertBoolean_2DTPP(tab)
    tab <- arrange_(tab, .dots = c(idVar, "temperature"))
    tab <- exportFct_sortCols(dat = tab, idVar = idVar, addCol = addCol, 
        intensityStr = intensityStr, fcStr = fcStr, normalizedData = normalizedData)
    if (addPlotColumns) {
        tab <- exportFct_addPlotColumns(tab = tab, path = outPath, 
            idVar = idVar, trRef = trRef)
    }
    tab <- exportFct_ensureUniqueColumns(tab)
    allCfgCols <- colnames(configTable)
    labelCols <- detectLabelColumnsInConfigTable(allColumns = allCfgCols)
    configTable[, labelCols] <- suppressWarnings(apply(configTable[, 
        labelCols], 2, as.numeric))
    wb <- createWorkbook()
    addWorksheet(wb, "Exp details")
    addWorksheet(wb, "pEC50")
    headerStyle <- createStyle(fgFill = "#DCE6F1", border = "Bottom", 
        textDecoration = "Bold")
    writeDataTable(wb, sheet = "Exp details", x = configTable, 
        startCol = 1, startRow = 1, rowNames = FALSE, colNames = TRUE, 
        headerStyle = headerStyle)
    writeDataTable(wb, sheet = "pEC50", x = tab, startCol = 1, 
        startRow = 1, rowNames = FALSE, colNames = TRUE, headerStyle = headerStyle)
    fc_cols <- grep(".*unmodified", colnames(tab))
    wb <- exportFct_colorCodeColumns(wb = wb, sheet = "pEC50", 
        cols = fc_cols, dat = tab)
    wb <- exportFct_addPlotLinks_2DTPP(wb = wb, sheet = "pEC50", 
        dat = tab)
    success <- exportFct_trySave(wb = wb, file = fileName)
    return(fileName)
}
<bytecode: 0x000000001bb27a40>
<environment: namespace:TPP>