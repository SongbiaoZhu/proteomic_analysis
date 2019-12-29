tpp2dCurveFit
function (configFile = NULL, data, nCores = 1, naStrs = NULL, 
    fcStr = NULL, idVar = NULL, nonZeroCols = NULL, r2Cutoff = 0.8, 
    fcCutoff = 1.5, slopeBounds = c(1, 50), fcTolerance = 0.1) 
{
    if (!missing(configFile)) {
        warning("`configFile` is deprecated.", call. = TRUE)
    }
    if (!missing(naStrs)) {
        warning("`naStrs` is deprecated.", call. = TRUE)
    }
    if (!missing(fcStr)) {
        warning("`fcStr` is deprecated.", call. = TRUE)
    }
    if (!missing(idVar)) {
        warning("`idVar` is deprecated.", call. = TRUE)
    }
    if (!missing(nonZeroCols)) {
        warning("`nonZeroCols` is deprecated.", call. = TRUE)
    }
    checkFunctionArgs(match.call(), c("data"))
    configTable <- attr(data, "configTable")
    importSettings <- attr(data, "importSettings")
    uniqueIdCol <- importSettings$uniqueIdCol
    if (is.null(uniqueIdCol)) {
        stop("attr(data, 'uniqueIdCol') must contain a field named 'uniqueIdCol'.")
    }
    else if (!is.character(uniqueIdCol)) {
        stop("attr(data, 'importSettings')$uniqueIdCol must be of class character.")
    }
    else {
        message("Looking for unique ID column: '", uniqueIdCol, 
            "'")
    }
    if (!uniqueIdCol %in% colnames(data)) {
        stop("Please specify an uniqueIdCol character string argument that represents a suffix of one of \n         the column names of your data!")
    }
    else if (length(data[[uniqueIdCol]]) != length(unique(data[[uniqueIdCol]]))) {
        stop("Please indicate an uniqueIdCol character string that matches a column with unique identifiers!")
    }
    nonZeroCols <- importSettings$nonZeroCols
    if (is.null(nonZeroCols)) {
        stop("attr(data, 'importSettings') must contain a field named 'nonZeroCols'.")
    }
    else if (!is.character(nonZeroCols)) {
        stop("attr(data, 'importSettings')$nonZeroCols must be of class character.")
    }
    else {
        message("Looking for nonZeroCols: '", nonZeroCols, "'")
    }
    if (!all(nonZeroCols %in% colnames(data))) {
        stop("The given QC columns (specified by attr(data, 'importSettings')$nonZeroCols) were not found in the column names of 'data'.")
    }
    finalFcPrefix <- obtain_fcStr_from_df_annotation(dat = data)
    message("Performing TPP-CCR dose response curve fitting and generating result table...")
    cfgIn <- convert_2D_cfgTable_to_CCR_cfgTable(configTable = configTable)
    datIn <- as.data.frame(data)
    CCRresult <- suppressMessages(analyzeTPPCCR(configTable = cfgIn, 
        data = datIn, nCores = nCores, resultPath = NULL, plotCurves = FALSE, 
        fcStr = finalFcPrefix, naStrs = c("NA", "n/d", "NaN", 
            "<NA>"), qualColName = "qupm", xlsxExport = FALSE, 
        idVar = uniqueIdCol, nonZeroCols = nonZeroCols, normalize = FALSE, 
        r2Cutoff = r2Cutoff, fcCutoff = fcCutoff, slopeBounds = slopeBounds, 
        fcTolerance = fcTolerance, ggplotTheme = NULL, verbose = FALSE))
    compound <- as.character(cfgIn$Experiment)
    allCols <- colnames(CCRresult)
    newCols <- sub(paste("*_", compound, sep = ""), "", allCols)
    colnames(CCRresult) <- newCols
    message("Done.")
    importSettings$r2Cutoff <- r2Cutoff
    importSettings$fcCutoff <- fcCutoff
    importSettings$slopeBounds <- slopeBounds
    importSettings$fcTolerance <- fcTolerance
    importSettings$uniqueIdCol <- "Protein_ID"
    attr(CCRresult, "importSettings") <- importSettings
    attr(CCRresult, "configTable") <- configTable
    return(CCRresult)
}
<bytecode: 0x0000000017990f00>
<environment: namespace:TPP>