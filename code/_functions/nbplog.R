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
#       x     = T by N matrix of data
#       kmax  = number of factors to estimate
#   
#   Outputs:
#   --------------------------------------------- 
#       IC2   =
#       index = 
# 
#---------------------------------------------------------------------------------------------------

nbplog <- function(x, kmax) {
    # Preliminaries
    IC2 <- matrix(0, kmax, 1)
    x <- scale(x)
    bigT <- nrow(x)
    bigN <- ncol(x)
    
    # SVD of Sums of Squares and Cross Products matrix
    xx <- x %*% t(x)
    fhat0 <- svd(xx)$u
    
    for (k in 1:kmax) {
        
        # Estimating factors and errors
        fhat <- as.matrix(fhat0[, 1:k]) * sqrt(bigT)
        lambda <- t(x) %*% fhat / bigT
        ehat <- x - fhat %*% t(lambda)
        
        # Panel information criteria for factor selection
        NT <- bigT * bigN
        NT1 <- bigT + bigN
        GCT <- min(bigT, bigN)
        
        # Computing the criterion
        CT2 <- log(GCT) * k * (NT1 / NT)
        
        # Compute value of information criterion
        Sigma <- mean() 
        
        
    }
}