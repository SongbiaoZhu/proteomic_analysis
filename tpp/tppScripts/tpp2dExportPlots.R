tpp2dExportPlots
function (plotList, resultPath, type = "none") 
{
    plotPath <- file.path(resultPath, "plots")
    if (!file.exists(plotPath)) {
        dir.create(file.path(resultPath, "plots"), recursive = TRUE)
    }
    invisible(lapply(names(plotList), function(pl) {
        if (!is.null(plotList[[pl]]) && !is.na(plotList[[pl]])) {
            if (type == "all" | type == "good") {
                savePl <- try(suppressMessages(ggsave(plotList[[pl]], 
                  filename = file.path(plotPath, paste(gsub("\\|", 
                    "_", gsub("(\\.)", "_", pl)), "2D_TPP", type, 
                    "plots.pdf", sep = "_")))))
                if (class(savePl) == "try-error") {
                  setwd(plotPath)
                  suppressMessages(ggsave(plotList[[pl]], filename = paste(gsub("\\|", 
                    "_", gsub("(\\.)", "_", pl)), "2D_TPP", type, 
                    "plots.pdf", sep = "_")))
                }
            } else if (type == "single") {
                if (!is.null(plotList[[pl]])) {
                  invisible(lapply(names(plotList[[pl]]), function(temp) {
                    savePl <- try(suppressMessages(ggsave(plotList[[pl]][[temp]], 
                      filename = file.path(plotPath, paste(gsub("\\|", 
                        "_", gsub("(\\.)", "_", pl)), gsub("(\\.)", 
                        "_", temp), "2D_TPP", type, "plots.pdf", 
                        sep = "_")))))
                    if (class(savePl) == "try-error") {
                      setwd(plotPath)
                      suppressMessages(ggsave(plotList[[pl]], 
                        filename = paste(gsub("\\|", "_", gsub("(\\.)", 
                          "_", pl)), gsub("(\\.)", "_", temp), 
                          "2D_TPP", type, "plots.pdf", sep = "_")))
                    }
                  }))
                }
            } else if (type == "spline") {
                savePl <- try(suppressMessages(ggsave(plotList[[pl]], 
                  filename = file.path(plotPath, paste(gsub("\\|", 
                    "_", gsub("(\\.)", "_", pl)), "2D_TPP", type, 
                    "plots.pdf", sep = "_")))))
                if (class(savePl) == "try-error") {
                  setwd(plotPath)
                  suppressMessages(ggsave(plotList[[pl]], filename = paste(gsub("\\|", 
                    "_", gsub("(\\.)", "_", pl)), "2D_TPP", type, 
                    "plots.pdf", sep = "_")))
                }
            } else {
                stop("Please specify a valid argument for 'type' ('all', 'good','single' or 'spline)!")
            }
        }
    }))
    return(NULL)
}
<bytecode: 0x000000001106a260>
<environment: namespace:TPP>