#### Quadratic Programming function for a subset of COGs ####
## M is a matrix with rows = COGs and columns = SCs
## E is a matrix with rows = COGs and columns = 1
## S is a vector with the names of the COGs subsetted 

subsetQP <- function(M, E, S){ ## M is the original matrix, E is the God frequency matrix e1, S is the COGs subset
  S <- data.frame(COG = S)
  E <- as.data.frame(E) %>% rownames_to_column(var = "COG") %>% 
    right_join(S, by="COG") %>% select(-COG) %>% as.matrix()
  M <- as.data.frame(M) %>% rownames_to_column(var = "COG") %>% 
    right_join(S, by="COG") %>% select(-COG) %>% as.matrix()
  output <- QP(M, E)
  return(output)
}