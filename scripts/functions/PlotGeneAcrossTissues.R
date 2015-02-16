# PlotGeneAcrossTissues.R

PlotGeneAcrossTissues <- function(dat, jtitle){
  m <- ggplot(dat, aes(x = time, y = exprs,
                       group = experiment, 
                       colour = experiment)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = "log2 expression")
  return(m)
}
