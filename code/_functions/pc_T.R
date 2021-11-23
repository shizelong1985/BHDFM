#---------------------------------------------------------------------------------------------------
#   This function estimates k latent common factors based on a panel of variables using the method 
#   of principal components. It is assumed that the panel can be described by a 
#   factor model of the form 
#           
#           X = F * Lambda' + e
#           
#   This function relies on the identifying normalization that F'F/T = I.
#   
#   Inputs:
#   --------------------------------------------- 
#       x       = T by N matrix of data
#       k       = number of factors to estimate
#   
#   Outputs:
#   --------------------------------------------- 
#       fhat    = T by k matrix of estimated factors
#       lambda  = N by k matrix of estimated factor loadings 
#       ehat    = residuals
# 
#---------------------------------------------------------------------------------------------------

pc_T <- function (x, k) {
    
    # Setting up preliminaries
    output  <- list()
    bigT    <- nrow(x)
    fhat    <- matrix(NA, nrow = bigT, ncol = k)
    
    # SVD of Sum of Squares and Cross Products matrix
    xx      <- x %*% t(x)   # T by T matrix
    fhat0   <- svd(xx)$u    # T by T matrix
    
    # Estimating factors, loadings, and errors
    fhat[, 1:k] <- fhat0[, 1:k] * sqrt(bigT)    # T by k matrix    
    lambda      <- t(x) %*% fhat / bigT         # N by k matrix
    ehat        <- x - fhat %*% t(lambda)       # T by N matrix
    
    output$fhat     <- fhat
    output$lambda   <- lambda
    output$ehat     <- ehat
    return(output)
}