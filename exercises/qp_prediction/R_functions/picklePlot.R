# picklePlot
# arguments: D=dataframe with COG frequency data;
#           X, Y = column names of datasets to compare;
#           PlotTitle = comparison name

picklePlot <- function(D, X, Y, PlotTitle){
  X <- enquo(X)
  Y <- enquo(Y)
  plot <- ggplot(D) +
    geom_point(aes(x=!!X, y=!!Y), alpha=0.3) +
    geom_abline(intercept=0, slope=1, lwd=0.5, lty=3, color="red") +
    xlim(0,1) + ylim(0,1) + theme_bw() + ggtitle(PlotTitle) +
    theme(plot.title=element_text(hjust=0.5))
  return(plot)
}