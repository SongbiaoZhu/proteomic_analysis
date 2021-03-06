rm(list = ls())
source("functions.R")
cre_folder("out_folder")

library(ggplot2)
theme_set(
  theme_gray() + theme(
    axis.line = element_line(size = 0.5),
    panel.background = element_rect(fill = NA, size = rel(20)),
    panel.grid.minor = element_line(colour = NA),
    axis.text = element_text(size = 10, colour = 'black'),
    axis.title = element_text(size = 14)
  )
)

if(file.exists("ipa_pathway.xlsx")){
  dat_ipapath <- openxlsx::read.xlsx("ipa_pathway.xlsx", sheet = 1)
  colnames(dat_ipapath) <- c("Pathway",
                             "pvalue",
                             "Ratio",
                             "zscore",
                             "abs.zscore",
                             "up",
                             "down")
  #Below are 3 fill color choice
  # rainbow color
  dat_ipapath = dat_ipapath[-c(6,19),]
  ggplot(data = dat_ipapath, aes(Pathway, pvalue)) +
    geom_bar(stat = "identity", aes(fill = zscore)) +
    labs(y = "-log10(p-value)", x = NULL) +
    scale_fill_gradientn(colours = rainbow(10)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip()
  # 2 color for negative and positive zscore
  ggplot(data = dat_ipapath, aes(Pathway, pvalue)) +
    geom_bar(stat = "identity", aes(fill = zscore < 0)) +
    labs(y = "-log10(p-value)", x = NULL) +
    scale_color_discrete(c('#1B01FF', '#FF6F00')) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip()
  # diveraging gradient color for negative and positive zscores
  ggplot(data = dat_ipapath, aes(Pathway, pvalue)) +
    geom_bar(stat = "identity", aes(fill = zscore)) +
    labs(y = "-log10(p-value)", x = NULL) +
    scale_fill_gradient2(
      low = "#1B01FF",
      mid = "white",
      high = "#FF6F00",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip()
  
  # Importantly, manually order the bar with biology
  ## set the levels in order we want
  ind <-
    rev(dat_ipapath[, 1][c(1:20)])
  dat_ipapath$Pathway <- factor(dat_ipapath$Pathway, levels = ind)
  
  ## Final plot
  my_y_title = expression(paste("Log"[10]*" "*italic("P"), " value"))
  # my_y_title <- expression(paste("-Log ", italic("P"), " value"))
  ggplot(data = dat_ipapath, aes(Pathway, pvalue)) +
    geom_bar(stat = "identity", aes(fill = zscore))+
    geom_text(aes(label=zscore), 
              position=position_dodge(width=0.9), 
              hjust=1,size = 2.5) +
    labs(y = my_y_title, x = NULL) +
    scale_fill_gradient2(
      name = 'Z-score',
      low = "#1B01FF",
      mid = "white",
      high = "#FF6F00",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    geom_vline(
      xintercept = c(4.5, 8.5, 13.5),
      linetype = "dashed",
      color = "black",
      size = 0.5
    ) + # coord_flip, the x y axis not change 
    annotate(
      geom = "text",
      x = c(18, 12, 7, 3),
      y = max(dat_ipapath$pvalue),
      label = c("Metabolism", "Immunity", "Apoptosis", "Motility"),
      size = 5,
      hjust= 1 #  hjust=1, is right aligned.
    ) 
  ggsave("out_folder/ipa_pathway_GPB.jpeg", 
         width = 170,height = 170,units = c('mm'))
  # save data
  save(dat_ipapath,ind, file = "output_step10.RData")
}else
  {stop("we can not find the file:ipa_pathway.xlsx")}







