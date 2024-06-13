#### Plot QP prediction outcome ####

# D is a dataframe with 4 columns:
# "GPSC" = strain classification, can be any classification doesn't have to necessarily be GPSC
# "vaccine" = vaccine type of the strain, either VT or NVT
# "SC_obs" = strain prevalences that were observed in the population
# "SC_pred" = strain prevalences that were predicted with the QP model

plotPred <- function(D, Title){
  SSE <- round(sum((D$SC_pred-D$SC_obs)^2), digits = 4)
  
  plot <- D %>% ggplot(aes(x=SC_pred, y=SC_obs, colour=vaccine)) + 
    geom_segment(aes(x=0,xend=0.14,y=0,yend=0.14), 
                 color="black",alpha=.7,lwd=0.5,lty=3) +
    theme(legend.position = "none") + theme_classic() + 
    geom_smooth(method='lm', color="#899DA4" ,
                formula=y~x, alpha=0.3, lwd=.6, 
                fullrange=T, linetype="blank", show.legend=F) +
    annotate(geom = "text", x=0.125, y =0.13, 
             label = "1:1 line", angle = 45, size = 3) + 
    geom_point(size=3, alpha = 0.9) + ggtitle(Title) + ## 
    scale_x_continuous("Predicted Prevalence")+
    scale_y_continuous("Observed Prevalence") +
    coord_fixed(ratio = 1, xlim=c(0,0.15), ylim=c(0,0.15)) + 
    scale_colour_manual(breaks=c("mixed", "NVT","VT"),
                        values = c("goldenrod", "#143c77","mediumorchid4")) + 
    theme(legend.position = "none", 
          axis.text = element_text(colour = "black"),) + 
    # geom_text_repel(aes(label = paste("SC", SC_global, sep = "-")), 
    #               data = outlier3B, size = 3.5) + 
    annotate("text", x= 0.015 , y = 0.15, size = 3, #col = "darkblue",
             label = paste("SSE =",SSE, sep = " "))
  return(plot)
}