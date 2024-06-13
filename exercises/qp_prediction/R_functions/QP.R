#### Quadratic Programming function ####
## M is a matrix with rows = COGs and columns = SCs
## E is a matrix with rows = COGs and columns = 1

QP <- function(M, E){ 
  #### check correct order ####
  x1 <- rownames(M)
  x2 <- rownames(E)
  test <- identical(x1,x2)
  if(test == FALSE) {
    "wrong order for QP/n"
    return(NULL)
  }
  
  rinv <- solve(chol(t(M) %*% M)) # M to be minimized in quad. function (Choleski decomp) [solve function here is giving the inverse matrix]
  C <- cbind(rep(1,ncol(M)), diag(ncol(M))) #Constraints to minimize the quad. function 
  b <- c(1,rep(0,ncol(M)))
  d <- t(E) %*% M  #Vector in the quadratic function to be minimized
  output <- solve.QP(Dmat = rinv, factorized = TRUE, dvec = d, 
                     Amat = C, bvec = b, meq = 1)$solution
  output <- round(output, digits = 5)
  return(output)
}