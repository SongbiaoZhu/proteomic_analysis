good_code
fTmp <- paste0(format(Sys.time(), "%Y-%m-%d"), "_results_2D_TPP.xlsx")
fileName <- file.path(outPath, fTmp)
message("Writing results to file: ", fileName)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Exp details")
headerStyle <- createStyle(border = "Bottom", 
                         textDecoration = "Bold")
writeDataTable(wb, sheet = "Exp details", x = config_tpp2d, 
             startCol = 1, startRow = 1, rowNames = FALSE, colNames = TRUE, 
             headerStyle = headerStyle)
saveWorkbook(wb, "writeDataTableExample.xlsx", overwrite = TRUE)

t1 <- Sys.time()
# run code here
timeDiff <- Sys.time() - t1
message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 
    2), " ", units(timeDiff), "\n")