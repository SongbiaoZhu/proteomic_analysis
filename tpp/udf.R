# user defined functions

creDir <- function(folder) {
  ifelse(!dir.exists(file.path(getwd(),
                               folder)), dir.create(file.path(getwd(), folder)), FALSE)
}

# create a function kegg_plot
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-Log10 P value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}

# convert first letter to upcase
firstup <- function(x) {
  x = tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# extract the calculate ratio
## "Abundance Ratio: (129, Sample) / (126, Control)" --to-- 'ratio129.126'
nameRatio <- function(x){
  name <- c()
  for (i in seq(length(x))) {
    y = paste0('ratio',paste(str_extract_all(x[i], "1\\d+")[[1]],collapse = '.'))
    name <- c(name,y)
  }
  return(name)
}
