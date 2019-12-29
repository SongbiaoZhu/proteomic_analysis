analyze2DTPP
function (configTable, data = NULL, resultPath = NULL, idVar = "gene_name", 
    fcStr = NULL, intensityStr = "signal_sum_", naStrs = c("NA", 
        "n/d", "NaN", "<NA>"), methods = "doseResponse", qualColName = "qupm", 
    compFc = TRUE, normalize = TRUE, addCol = NULL, nCores = 1, 
    nonZeroCols = "qssm", fcTolerance = 0.1, r2Cutoff = 0.8, 
    fcCutoff = 1.5, slopeBounds = c(1, 50), fractAbund = FALSE, 
    xlsxExport = TRUE, plotAll = FALSE, plotAllR2 = FALSE, plotSingle = FALSE, 
    trRef = NULL, refFcStr = "norm_rel_fc_", addInfo = FALSE, 
    createReport = "none", paletteName = "Spectral", configFile) 
{
    message("This is TPP version ", packageVersion("TPP"), ".")
    if (!missing(configFile)) {
        warning("`configFile` is deprecated. Use 'configTable' instead.", 
            call. = FALSE)
        configTable <- configFile
    }
    datIn <- tpp2dImport(configTable = configTable, data = data, 
        idVar = idVar, addCol = addCol, intensityStr = intensityStr, 
        qualColName = qualColName, nonZeroCols = nonZeroCols, 
        fcStr = fcStr)
    confgFields <- suppressMessages(importCheckConfigTable(infoTable = configTable, 
        type = "2D"))
    files <- suppressWarnings(confgFields$Path)
    outDirList <- importFct_makeOutputDirs(outDir = resultPath, 
        fNames = files)
    flagDoWrite <- outDirList$doWrite
    resultPath <- outDirList$outDir
    if (!flagDoWrite) 
        xlsxExport = plotAll = plotAllR2 = plotSingle <- FALSE
    if (compFc) {
        fcStr <- "rel_fc_protein_"
        datIn <- tpp2dComputeFoldChanges(data = datIn, newFcStr = fcStr)
    }
    if (normalize) {
        NormData2d <- tpp2dNormalize(data = datIn)
    }
    else {
        NormData2d <- datIn
    }
    fcStrUpdated <- attr(NormData2d, "importSettings")$fcStrNorm
    if (length(which(is.na(NormData2d[[qualColName[1]]]))) != 
        0) {
        NormData2d <- NormData2d[-which(is.na(NormData2d[[qualColName[1]]])), 
            ]
    }
    if (length(which(duplicated(NormData2d$unique_ID))) != 0) {
        message("There are duplicated proteins in your experimental conditions! These are filtered out to run this analysis!\nPlease check your data quality and consider pre-filtering!")
        NormData2d <- NormData2d[-which(duplicated(NormData2d$unique_ID)), 
            ]
    }
    if (fractAbund) {
        NormData2d <- tpp2dCalcFractAbundance(data = NormData2d)
    }
    if ("doseResponse" %in% methods) {
        analysisResults <- tpp2dCurveFit(data = NormData2d, nCores = nCores, 
            r2Cutoff = r2Cutoff, fcCutoff = fcCutoff, slopeBounds = slopeBounds, 
            fcTolerance = fcTolerance)
        if (plotAll) {
            plotList <- tpp2dCreateDRplots(data = analysisResults, 
                type = "all", verbose = TRUE, paletteName = paletteName)
            tpp2dExportPlots(plotList = plotList, resultPath = resultPath, 
                type = "all")
        }
        if (plotAllR2) {
            plotGoodList <- tpp2dCreateDRplots(data = analysisResults, 
                type = "good", verbose = TRUE, paletteName = paletteName)
            tpp2dExportPlots(plotList = plotGoodList, resultPath = resultPath, 
                type = "good")
        }
        if (plotSingle) {
            plotSingleList <- tpp2dCreateDRplots(data = analysisResults, 
                type = "single", verbose = TRUE)
            tpp2dExportPlots(plotList = plotSingleList, resultPath = resultPath, 
                type = "single")
        }
    }
    else {
        analysisResults <- NormData2d
    }
    if (("splineFit" %in% methods) && !is.null(trRef)) {
        analysisResults <- tpp2dSplineFitAndTest(data = analysisResults, 
            dataRef = trRef, refIDVar = "Protein_ID", refFcStr = "norm_rel_fc_", 
            resultPath = resultPath, doPlot = TRUE, verbose = TRUE, 
            nCores = nCores)
    }
    else if ("splineFit" %in% methods) {
        message("The spline fit and corresponding f-Test could not be performed, as no TPP-TR reference dataset was specified!\nPlease check the file path you have specified for trRef!")
    }
    if (addInfo) {
        analysisResults <- tpp2dAddAdditionalInfo(data = analysisResults, 
            idVar = idVar)
    }
    if (!is.null(trRef) && file.exists(trRef)) {
        analysisResults <- tpp2dMerge2dRef(data = analysisResults, 
            trRef = trRef, idVar = idVar)
    }
    if (!is.null(resultPath) & xlsxExport) {
        addPlotColumns <- any(c(plotAll, plotAllR2, plotSingle, 
            !is.null(trRef)))
        tpp2dExport(tab = analysisResults, outPath = resultPath, 
            addCol = addCol, addPlotColumns = addPlotColumns, 
            trRef = trRef)
    }
    if ((createReport != "none") && !is.null(resultPath)) {
        tpp2dCreateReport(resultPath = resultPath, configFile = configTable, 
            normalize = normalize, configTable = configTable, 
            data = analysisResults, idVar = "Protein_ID", fcStr = fcStr, 
            fcStrUpdated = fcStrUpdated, documentType = createReport, 
            intensityStr = intensityStr, addCol = addCol, fcTolerance = fcTolerance, 
            r2Cutoff = r2Cutoff, fcCutoff = fcCutoff, slopeBounds = slopeBounds, 
            trRef = trRef)
    }
    return(analysisResults)
}
<bytecode: 0x0000000017ca7e08>
<environment: namespace:TPP>