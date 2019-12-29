tpptrAnalyzeMeltingCurves
function (data, pValMethod = "robustZ", pValFilter = list(minR2 = 0.8, 
    maxPlateau = 0.3), pValParams = list(binWidth = 300)) 
{
    message("Starting melting curve analysis.")
    expInfo <- sapply(data, annotation)
    dataSplit <- retrieveDataFromESets_TR(data)
    curveParDF <- dataSplit$curveParDF
    compDF <- createComparisonTable(infoTable = expInfo)
    if (!is.null(compDF)) {
        pValCols <- pValFctPerformComparisons(curveParsAllExp = curveParDF, 
            method = pValMethod, controlFilter = pValFilter, 
            controlpVal = pValParams, comparisons = compDF)
        qualCheckCols <- checkResultCols_tpptr(pValDF = pValCols, 
            curveParDF = curveParDF, comparisons = compDF)
    }
    else {
        pValCols <- qualCheckCols <- NULL
    }
    resultTable <- mergeOutputTables_TR(dataList = dataSplit, 
        pValDF = pValCols, qualCheckDF = qualCheckCols)
    return(resultTable)
}
<bytecode: 0x000000001458b9f0>
<environment: namespace:TPP>