analyzeTPPTR
function (configTable, data = NULL, resultPath = NULL, methods = c("meltcurvefit", 
    "splinefit"), idVar = "gene_name", fcStr = "rel_fc_", ciStr = NULL, 
    naStrs = c("NA", "n/d", "NaN", "<NA>"), qualColName = "qupm", 
    normalize = TRUE, normReqs = tpptrDefaultNormReqs(), ggplotTheme = tppDefaultTheme(), 
    nCores = "max", startPars = c(Pl = 0, a = 550, b = 10), splineDF = c(3:7), 
    maxAttempts = 500, plotCurves = TRUE, fixedReference = NULL, 
    pValMethod = "robustZ", pValFilter = list(minR2 = 0.8, maxPlateau = 0.3), 
    pValParams = list(binWidth = 300), verbose = FALSE, xlsxExport = TRUE) 
{
    meltcurve_plot = Protein_ID <- NULL
    message("This is TPP version ", packageVersion("TPP"), ".")
    trData <- tpptrImport(configTable = configTable, data = data, 
        idVar = idVar, fcStr = fcStr, naStrs = naStrs, qualColName = qualColName)
    if (!is.null(ciStr)) {
        options(TPPTR_CI = TRUE)
        trDataCI <- tpptrImport(configTable = configTable, data = data, 
            idVar = idVar, fcStr = ciStr, naStrs = naStrs, qualColName = qualColName)
        stopifnot(all(sapply(names(trData), function(i) {
            all(dim(trData[[i]]) == dim(trDataCI[[i]]))
        })))
    }
    else {
        options(TPPTR_CI = FALSE)
        trDataCI <- NULL
    }
    expInfo <- sapply(trData, annotation)
    expNames <- names(trData)
    expNum <- length(expNames)
    expConds <- sapply(trData, function(s) s@annotation[["condition"]])
    expComps <- createComparisonTable(infoTable = expInfo)
    confgFields <- suppressMessages(importCheckConfigTable(infoTable = configTable, 
        type = "TR"))
    files <- confgFields$files
    outDirList <- importFct_makeOutputDirs(outDir = resultPath, 
        fNames = files)
    flagDoWrite <- outDirList$doWrite
    pathDataObj <- outDirList$pathDataObj
    pathExcel = resultPath <- outDirList$outDir
    if (!flagDoWrite) 
        plotCurves <- FALSE
    if (flagDoWrite) {
        save(list = c("trData", "trDataCI"), file = file.path(pathDataObj, 
            "importedData.RData"))
    }
    if (normalize) {
        normResults <- tpptrNormalize(data = trData, normReqs = normReqs, 
            qcPlotTheme = ggplotTheme, qcPlotPath = NULL, fixedReference = fixedReference)
        trDataNormalized <- normResults[["normData"]]
        if (flagDoWrite) {
            save(list = c("trDataNormalized"), file = file.path(pathDataObj, 
                "normalizedData.RData"))
        }
    }
    else {
        trDataNormalized <- trData
    }
    allIDs <- lapply(trDataNormalized, featureNames) %>% unlist %>% 
        unname %>% unique %>% sort
    resultTable <- data.frame(Protein_ID = allIDs, stringsAsFactors = FALSE)
    if (any(methods == "splinefit")) {
        fitData <- trDataNormalized %>% tpptrTidyUpESets(., returnType = "exprs")
        if (!any(methods == "meltcurvefit")) {
            annotData <- trDataNormalized %>% tpptrTidyUpESets(., 
                returnType = "featureData")
        }
        else {
            annotData <- NULL
        }
        splineFitResultTable <- tpptrSplineFitAndTest(data = fitData, 
            factorsH1 = "condition", resultPath = resultPath, 
            doPlot = plotCurves, splineDF = splineDF, nCores = nCores, 
            additionalCols = annotData) %>% mutate(Protein_ID = as.character(Protein_ID))
        if (flagDoWrite) {
            save(list = c("splineFitResultTable"), file = file.path(pathDataObj, 
                "trResultsSplineFit.RData"))
        }
        resultTable <- left_join(resultTable, splineFitResultTable)
    }
    if (any(methods == "meltcurvefit")) {
        trDataFitted <- tpptrCurveFit(data = trDataNormalized, 
            dataCI = trDataCI, resultPath = resultPath, ggplotTheme = ggplotTheme, 
            doPlot = plotCurves, startPars = startPars, maxAttempts = maxAttempts, 
            nCores = nCores, verbose = verbose)
        if (flagDoWrite) {
            save(list = c("trDataFitted"), file = file.path(pathDataObj, 
                "fittedData.RData"))
        }
        meltCurveResultTable <- tpptrAnalyzeMeltingCurves(data = trDataFitted, 
            pValMethod = pValMethod, pValFilter = pValFilter, 
            pValParams = pValParams) %>% rename(meltcurve_plot = plot) %>% 
            mutate(meltcurve_plot = as.character(meltcurve_plot)) %>% 
            mutate(Protein_ID = as.character(Protein_ID))
        if (flagDoWrite) {
            save(list = c("meltCurveResultTable"), file = file.path(pathDataObj, 
                "trResultsMeltCurveFit.RData"))
        }
        resultTable <- suppressMessages(left_join(resultTable, 
            meltCurveResultTable))
    }
    if (flagDoWrite) {
        save(list = c("resultTable"), file = file.path(pathDataObj, 
            "results_TPP_TR.RData"))
    }
    if (xlsxExport) {
        if (flagDoWrite) {
            if (expNum > 1) {
                compNums <- assignCompNumber_to_expName(compDF = expComps, 
                  expNames = expNames)
                expColors <- plotColors(expConditions = expConds, 
                  comparisonNums = compNums)
            }
            else {
                expColors <- NULL
            }
            tppExport(tab = resultTable, file = file.path(pathExcel, 
                "results_TPP_TR.xlsx"), expColors = expColors, 
                expNames = expNames)
        }
        else {
            message("Cannot produce xlsx output because no result path is specified.")
        }
    }
    if (flagDoWrite) {
        pdf(file = file.path(resultPath, "QCplots.pdf"), width = 8, 
            height = 9)
        message("Creating QC plots to visualize median curve fits...")
        if (normalize) {
            qcPlotMedianFit <- normResults[["qcPlotObj"]]
            suppressWarnings(print(qcPlotMedianFit))
        }
        message("done.\n")
        message("Creating QC plots to visualize normalization effects...")
        qcPlotCorrelateGroupsRaw <- tppQCPlotsCorrelateExperiments(tppData = trData, 
            annotStr = "Non-normalized Fold Changes", ggplotTheme = ggplotTheme)
        if (normalize) {
            qcPlotCorrelateGroupsNorm <- tppQCPlotsCorrelateExperiments(tppData = trDataNormalized, 
                annotStr = "Normalized Fold Changes", ggplotTheme = ggplotTheme)
        }
        for (pn in names(qcPlotCorrelateGroupsRaw)) {
            suppressWarnings(print(qcPlotCorrelateGroupsRaw[[pn]]))
            if (normalize) 
                suppressWarnings(print(qcPlotCorrelateGroupsNorm[[pn]]))
        }
        message("done.\n")
        if (any(methods == "meltcurvefit")) {
            tpptrQCplots(resultTab = resultTable, expNames = expNames, 
                expConditions = expConds, compDF = expComps, 
                minR2 = pValFilter$minR2, ggplotTheme = ggplotTheme)
        }
        dev.off()
    }
    invisible(resultTable)
}
<bytecode: 0x00000000147d0b00>
<environment: namespace:TPP>