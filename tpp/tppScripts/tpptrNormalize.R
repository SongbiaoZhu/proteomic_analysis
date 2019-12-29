tpptrNormalize
function (data, normReqs = tpptrDefaultNormReqs(), qcPlotTheme = tppDefaultTheme(), 
    qcPlotPath = NULL, startPars = c(Pl = 0, a = 550, b = 10), 
    maxAttempts = 1, fixedReference = NULL) 
{
    infoNormP <- filterTables(data = data, normReqs = normReqs)
    normP <- infoNormP[["protein_IDs"]]
    listNormP <- sapply(data, function(x) {
        exprSubset(x, subset = normP)
    }, simplify = FALSE)
    message("-----------------------------------")
    listNormFit <- computeNormFactors(data = listNormP, startPars = startPars, 
        maxAttempts = maxAttempts, fixedReference = fixedReference)
    meltCurveModels <- listNormFit[["models"]]
    fcMedians <- listNormFit[["medians"]]
    tempVals <- listNormFit[["tempVals"]]
    r2Vec <- listNormFit[["rSquared"]]
    nNormP <- length(normP)
    qcPlotsMedianFits <- plotNormCurves(modelList = meltCurveModels, 
        xMat = tempVals, fcMat = fcMedians, r2Vec = r2Vec, nNormP = nNormP, 
        plotTheme = qcPlotTheme)
    if (!is.null(qcPlotPath)) {
        pdf(file = file.path(qcPlotPath, "QCplots_median_fits.pdf"), 
            width = 8, height = 9)
        print(qcPlotsMedianFits)
        dev.off()
    }
    message("-----------------------------------")
    message("Normalizing all proteins in all experiments.")
    dfCoeffs = listNormFit[["corrFactors"]]
    normData <- sapply(names(data), simplify = FALSE, function(n) {
        applyCoeffs(data = data[[n]], coeffs = dfCoeffs[, n])
    })
    message("Normalization successfully completed!\n")
    return(list(normData = normData, qcPlotObj = qcPlotsMedianFits))
}
<bytecode: 0x0000000013efd178>
<environment: namespace:TPP>