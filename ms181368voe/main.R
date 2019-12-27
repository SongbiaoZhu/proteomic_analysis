rm(list=ls())

# edit the step1.R, change the dataset name
scripts <- c("step1_data_clean.R", "step2_qc_identification.R",
             "step3_qc_quantification.R", "step4_filter_deps.R",
             "step5_export_deps_tables.R", "step6_volcano_text.R",
             "step7_gopathway.R", "step8_gsea.R",
             "step9_TFmap.R")
sapply(scripts, source) 
