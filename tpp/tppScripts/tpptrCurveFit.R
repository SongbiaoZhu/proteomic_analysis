tpptrCurveFit
function (data, dataCI = NULL, resultPath = NULL, ggplotTheme = tppDefaultTheme(), 
    doPlot = TRUE, startPars = c(Pl = 0, a = 550, b = 10), maxAttempts = 500, 
    nCores = "max", verbose = FALSE) 
{
    p <- NULL
    doPlot <- doPlot && !is.null(resultPath)
    useCI <- getOption("TPPTR_CI")
    if (is.null(useCI)) 
        useCI <- FALSE
    expInfo <- sapply(data, annotation)
    expNames <- expInfo["name", ]
    grConditions <- expInfo["condition", ]
    compDF <- createComparisonTable(infoTable = expInfo)
    protIDs <- unique(unlist(lapply(data, featureNames)))
    if (doPlot) {
        plotDir <- "Melting_Curves"
        if (!file.exists(file.path(resultPath, plotDir))) {
            dir.create(file.path(resultPath, plotDir), recursive = TRUE)
        }
        fNames <- paste0("meltCurve_", gsub("([^[:alnum:]])", 
            "_", protIDs))
        maxLen <- 255 - nchar(".pdf") - nchar("_truncated") - 
            nchar(as.character(length(fNames)))
        tooLong <- nchar(fNames) > maxLen
        cropSuffix <- paste0("_truncated", 1:sum(tooLong))
        fNames <- sapply(fNames, function(fTmp) {
            fNew <- substr(fTmp, 1, min(maxLen, nchar(fTmp)))
        }, simplify = TRUE, USE.NAMES = FALSE)
        fNames[tooLong] <- paste0(fNames[tooLong], cropSuffix)
        fNames <- paste0(fNames, ".pdf")
        plotPathsFull <- file.path(plotDir, fNames)
    }
    else {
        plotPathsFull <- rep("", length(protIDs))
    }
    names(plotPathsFull) <- protIDs
    xMat <- t(sapply(data, function(set) {
        set$temperature
    }))
    yMat <- data.frame(matrix(nrow = 0, ncol = (ncol(data[[1]]) + 
        2)))
    for (g in expNames) {
        yTmp <- Biobase::exprs(data[[g]])
        if (nrow(yTmp) == 0) {
            yTmp <- matrix(NA_real_, nrow = length(protIDs), 
                ncol = ncol(data[[g]]), dimnames = list(protIDs, 
                  colnames(yTmp)))
        }
        colnames(yTmp) <- 1:ncol(yTmp)
        yMat <- rbind(yMat, data.frame(expName = g, protID = rownames(yTmp), 
            FC = yTmp, stringsAsFactors = FALSE))
    }
    if (useCI) {
        ciMat <- data.frame(matrix(nrow = 0, ncol = (ncol(data[[1]]) + 
            2)))
        for (g in expNames) {
            ciTmp <- Biobase::exprs(dataCI[[g]])
            if (nrow(ciTmp) == 0) {
                ciTmp <- matrix(NA_real_, nrow = length(protIDs), 
                  ncol = ncol(dataCI[[g]]), dimnames = list(protIDs, 
                    colnames(ciTmp)))
            }
            colnames(ciTmp) <- 1:ncol(ciTmp)
            ciMat <- rbind(ciMat, data.frame(expName = g, protID = rownames(ciTmp), 
                CI = ciTmp, stringsAsFactors = FALSE))
        }
    }
    addLegend <- checkIfLegendPossible()
    message("Fitting melting curves to ", length(protIDs), " proteins.")
    nCores <- checkCPUs(cpus = nCores)
    t1 <- Sys.time()
    if (nCores > 1) {
        doParallel::registerDoParallel(cores = nCores)
        dfCurvePars <- foreach(p = protIDs, .combine = rbind, 
            .inorder = FALSE, .verbose = FALSE, .packages = c("ggplot2", 
                "gridExtra")) %dopar% {
            tpptrHelperFitandPlot(p, yMat, xMat, ciMat, startPars, 
                maxAttempts, expNames, verbose, ggplotTheme, 
                grConditions, compDF, addLegend, resultPath, 
                plotPathsFull, useCI, doPlot)
        }
        stopImplicitCluster()
    }
    else {
        dfCurvePars <- lapply(protIDs, function(p) {
            tpptrHelperFitandPlot(p, yMat, xMat, ciMat, startPars, 
                maxAttempts, expNames, verbose, ggplotTheme, 
                grConditions, compDF, addLegend, resultPath, 
                plotPathsFull, useCI, doPlot)
        })
        dfCurvePars = do.call("rbind", dfCurvePars)
    }
    timeDiff <- Sys.time() - t1
    message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 
        2), " ", units(timeDiff), "\n")
    gc(verbose = FALSE)
    message("Melting curves fitted sucessfully!")
    dfCurvePars$sufficient_data_for_fit <- as.logical(dfCurvePars$sufficient_data_for_fit)
    dfCurvePars$model_converged <- as.logical(dfCurvePars$model_converged)
    conv <- dfCurvePars$model_converged
    expr <- dfCurvePars$sufficient_data_for_fit
    message(sum(conv, na.rm = TRUE), " out of ", sum(expr), " models with sufficient data points converged (", 
        round(sum(conv, na.rm = TRUE)/sum(expr) * 100, 2), " %).\n")
    data <- storeMeltCurveParams(data = data, params = dfCurvePars)
    return(data)
}
<bytecode: 0x0000000013e1ad98>
<environment: namespace:TPP>