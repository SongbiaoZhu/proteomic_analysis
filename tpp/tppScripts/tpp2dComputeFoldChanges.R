tpp2dComputeFoldChanges
function (configTable = NULL, data, intensityStr = NULL, fcStr = NULL, 
    newFcStr = "rel_fc_") 
{
    if (!missing(configTable)) {
        warning("`configTable` is deprecated.", call. = TRUE)
    }
    if (!missing(intensityStr)) {
        warning("`intensityStr` is deprecated.", call. = TRUE)
    }
    if (!missing(fcStr)) {
        warning("`fcStr` is deprecated.", call. = TRUE)
    }
    checkFunctionArgs(match.call(), c("data"))
    configTable <- attr(data, "configTable")
    importSettings <- attr(data, "importSettings")
    intensityStr <- importSettings$intensityStr
    if (is.null(intensityStr)) {
        stop("attr(data, 'importSettings') must contain a field named 'intensityStr'.")
    }
    else if (!is.character(intensityStr)) {
        stop("attr(data, 'importSettings')$intensityStr must be of class character.")
    }
    else {
        message("Looking for intensity column prefix: '", intensityStr, 
            "'")
    }
    intensity.col.ids <- grep(intensityStr, colnames(data))
    if (length(intensity.col.ids) == 0) {
        stop("The given prefix for intensity columns (specified by attr(data, 'importSettings')$intensityStr) was not found in the column names of 'data'.")
    }
    message("Computing fold changes...")
    fc.cols <- sapply(intensity.col.ids, function(sc) {
        return(paste(newFcStr, sub(intensityStr, "", colnames(data)[sc]), 
            sep = ""))
    })
    ref.col.id <- grep(paste(intensityStr, "0$", sep = ""), colnames(data))
    if (length(ref.col.id) == 0) {
        ref.col.id <- grep(paste(intensityStr, "0.0$", sep = ""), 
            colnames(data))
        if (length(ref.col.id) == 0) {
            ref.col.id <- grep(paste(intensityStr, "0.00$", sep = ""), 
                colnames(data))
            if (length(ref.col.id) == 0) {
                ref.col.id <- grep(paste(intensityStr, "0.000$", 
                  sep = ""), colnames(data))
            }
        }
    }
    ref.vals <- as.numeric(as.character(data[, ref.col.id]))
    fc.m <- data.frame(matrix(as.numeric(as.character(unlist(data[, 
        intensity.col.ids]))), nrow = nrow(data[, intensity.col.ids]))/ref.vals)
    names(fc.m) <- fc.cols
    data[fc.cols] <- fc.m
    importSettings$fcStr <- newFcStr
    attr(data, "importSettings") <- importSettings
    message("Done.")
    return(data)
}
<bytecode: 0x0000000017358cc0>
<environment: namespace:TPP>