tpp2dCalcFractAbundance
function (configTable = NULL, data, intensityStr = NULL, idVar = NULL) 
{
    if (!missing(configTable)) {
        warning("`configTable` is deprecated.", call. = TRUE)
    }
    if (!missing(intensityStr)) {
        warning("`intensityStr` is deprecated.", call. = TRUE)
    }
    if (!missing(idVar)) {
        warning("`idVar` is deprecated.", call. = TRUE)
    }
    checkFunctionArgs(match.call(), c("data"))
    message("Calculating fractional abundance and DMSO1 vs. DMSO2 ratio...")
    configTable <- attr(data, "configTable")
    importSettings <- attr(data, "importSettings")
    intensityStr <- importSettings$intensityStr
    idVar <- importSettings$proteinIdCol
    compound <- as.character(configTable$Compound[1])
    all.temps <- configTable$Temperature
    temps.num <- length(all.temps)
    temp.pairs <- split(all.temps, as.factor(sort(rank(all.temps)%%(temps.num/2))))
    intensity.cols <- grep(intensityStr, colnames(data))
    subset.list <- lapply(unique(data[[idVar]]), function(prot) {
        CCR.subset <- data[which(data[[idVar]] == prot), ]
        if (nrow(CCR.subset) == temps.num) {
            temp.ratios <- sapply(temp.pairs, function(t.pair) {
                row1 <- grep(as.character(t.pair[1]), CCR.subset[["temperature"]])
                row2 <- grep(as.character(t.pair[2]), CCR.subset[["temperature"]])
                dmso.ratio <- as.numeric(as.character(CCR.subset[row1, 
                  intensity.cols[length(intensity.cols)]]))/as.numeric(as.character(CCR.subset[row2, 
                  intensity.cols[length(intensity.cols)]]))
                f.abund1 <- 100 * sum(as.numeric(sapply(CCR.subset[row1, 
                  intensity.cols], as.character)))/sum(as.numeric(sapply(CCR.subset[row1, 
                  intensity.cols], as.character)), as.numeric(sapply(CCR.subset[row2, 
                  intensity.cols], as.character)))
                f.abund2 <- 100 * sum(as.numeric(sapply(CCR.subset[row2, 
                  intensity.cols], as.character)))/sum(as.numeric(sapply(CCR.subset[row1, 
                  intensity.cols], as.character)), as.numeric(sapply(CCR.subset[row2, 
                  intensity.cols], as.character)))
                return(list(dmso.ratio, f.abund1, f.abund2))
            })
            CCR.subset$dmso1_vs_dmso2 <- unlist(rep(temp.ratios[1, 
                ], each = 2))
            CCR.subset$total_sumionarea_fraction <- as.vector(rbind(unlist(temp.ratios[2, 
                ]), unlist(temp.ratios[3, ])))
        }
        else {
            CCR.subset$dmso1_vs_dmso2 <- rep(NA, nrow(CCR.subset))
            CCR.subset$total_sumionarea_fraction <- rep(NA, nrow(CCR.subset))
        }
        return(CCR.subset)
    })
    out <- do.call(rbind, subset.list)
    attr(out, "importSettings") <- importSettings
    attr(out, "configTable") <- configTable
    return(out)
}
<bytecode: 0x0000000018fd33f8>
<environment: namespace:TPP>