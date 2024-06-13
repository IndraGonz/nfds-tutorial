#### Plot to get SSE number from prediction ####

# D is a dataframe with 4 columns:
# "GPSC" = strain classification, can be any classification doesn't have to necessarily be GPSC
# "vaccine" = vaccine type of the strain, either VT or NVT
# "SC_obs" = strain prevalences that were observed in the population
# "SC_pred" = strain prevalences that were predicted with the QP model

SSEpred <- function(D){
  SSE <- round(sum((D$SC_pred-D$SC_obs)^2), digits = 4)
  return(SSE)
}