#---------------------------------------------------------------------------------------------------
#   
#   Gets principal component analysis (PCA) estimates for sub-block-specific factors, block-specific
#   factors, and common factors of a hierarchical dynamic factor model of the form:
#   
#           X(t) = Lambda_H * H(t) + e
#           H(t) = Lambda_G * G(t) + e
#           G(t) = Lambda_F * F(t) + e
#           F(t) = Psi * F(t-1) + e
#   
#   If there is a sub-block structure in the data, use PCA on the sub-block data to get estimates of
#   the sub-block factors, Then, use PCA on the estimated sub-block factors to get estimates of the 
#   block factors. Finally, use PCA on the estimated block factors to get estimates of the common 
#   factors. then use the sub-block factor estimate. If there is only a block structure in the data, 
#   use PCA on the block data to get estimates of the block factors. Then, use PCA on the estimated
#   block factors to get estimates of the common factors.
#           
#   
#   Inputs:
#   --------------------------------------------- 
#       X    = list of block/sub-block data.
#       Bsub = vector consisting of the number of sub-blocks in each block
#       K_F  = number of common factors
#       K_G  = vector consisting of the number of block-specific factors in each block
#       K_H  = vector consisting of the number of sub-block-specific factors in each sub-block
#   
#   Outputs:
#   --------------------------------------------- 
#       H_PC = PC estimates of sub-block factors for each sub-block
#       G_PC = PC estimates of block factors for each block
#       F_PC = PC estimates of common factors
# 
#---------------------------------------------------------------------------------------------------

do_PC_level4 <- function(X, Bsub, K_F, K_G, K_H) {
    
    # Setting up preliminaries
    output <- list()
    B <- length(Bsub)
    bigT <- nrow(X[[1]][[1]])
    bigGmat <- matrix(NA, nrow = bigT, ncol = 0)
    bigHmat <- matrix(NA, nrow = bigT, ncol = 0)
    G_PC <- vector("list", B)
    H_PC <- vector("list", sum(Bsub != 0))
    F_PC <- vector("double", bigT)
    
    for (b in 1:B) {
        # Xb is all data from block b
        Xb <- X[[b]]
        if (is.list(Xb)) {
            Xb <- matrix(unlist(Xb), nrow = bigT)
        }
        
        for (s in 1:Bsub[b]) {
            if (Bsub[b] > 0) { # If Bsub[b] > 0, then there is sub-block structure
                # Zbs is all data from sub-block s of block b
                Zbs <- X[[b]][[s]]
                H_PC[[b]][[s]] <- pc_T(scale(Zbs), K_H[[b]][[s]])$fhat
                bigHmat <- cbind(bigHmat, H_PC[[b]][[s]])
                
                G_PC[[b]] <- pc_T(scale(matrix(unlist(H_PC[[b]]), nrow = bigT)), K_G[b])$fhat
            } else {
                G_PC[[b]] <- pc_T(scale(Xb), K_G[b])$fhat
            }
        }
        bigGmat <- cbind(bigGmat, G_PC[[b]])
    } 
    
    F_PC <- pc_T(bigGmat, K_F)$fhat
    
    output$H_PC = H_PC
    output$G_PC = G_PC
    output$F_PC = F_PC
    
    return(output)
}
