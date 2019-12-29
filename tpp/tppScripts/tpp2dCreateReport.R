tpp2dCreateReport
function (data = NULL, configFile = NULL, resultPath = NULL, 
    documentType = "html_document", configTable = NULL, normalize = TRUE, 
    methods = c(""), idVar = "gene_name", fcStr = "rel_fc_", 
    fcStrUpdated = "norm_rel_fc_", intensityStr = "signal_sum_", 
    addCol = NULL, fcTolerance = NA, r2Cutoff = NA, fcCutoff = NA, 
    slopeBounds = c(NA, NA), fTest = FALSE, trRef = "none") 
{
    resultTable <- data
    if (!is.null(resultPath)) {
        if (!file.exists(resultPath)) {
            dir.create(resultPath, recursive = TRUE)
        }
        message("Creating report...\n")
        inFile <- file.path(system.file(package = "TPP"), "tpp2d_report.Rmd")
        if (documentType == "html_document") {
            outFile <- file.path(resultPath, paste("report_", 
                basename(resultPath), ".html", sep = ""))
        }
        else if (documentType == "pdf_document") {
            outFile <- file.path(resultPath, paste("report_", 
                basename(resultPath), ".pdf", sep = ""))
        }
        else {
            stop("Please specify a valid documentType!")
        }
        render(input = inFile, output_file = outFile, envir = environment(), 
            intermediates_dir = tempdir(), output_format = documentType)
    }
    else {
        warning("Report could not be produced because resultPath was not specified.")
    }
}
<bytecode: 0x0000000017479470>
<environment: namespace:TPP>