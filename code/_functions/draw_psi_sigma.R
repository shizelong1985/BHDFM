#---------------------------------------------------------------------------------------------------
#   
#   Inputs:
#   --------------------------------------------- 
#   
#   
#   
#   Outputs:
#   --------------------------------------------- 
#       psi   = vector of length K*q consisting of autoregressive coefficients      
#       sigma = 
#   
#---------------------------------------------------------------------------------------------------

draw_psi_sigma <- function (Ft, psi_F, sig2_F, prior, params, var_fix, ...) {

    # Preliminaries
    bigT <- nrow(Ft)
    K_F <- ncol(Ft)
    l_F <- getElement(params, "l_F")
    q_F <- getElement(params, "q_F")
    max_q_F <- max(q_F)
    max_l_F <- max(l_F)
    output <- list()
    
    PSI <- psi_F  # K_F by 1
    SIG <- sig2_F # K_F by 1
    
    for (k in 1:K_F) {
        q_Fk <- q_F[k]
        xpsi_F1 <- matrix(0, nrow = max_q_F, 1)
        
        if (q_Fk > 0) {
            ystar <- as.matrix(Ft[, k]) # bigT by K_F
            xstar <- NA
            
            for (i in 1:q_Fk) {
                xstar <- cbind(xstar, lagn(Ft[, k], i)) # bigT by K_F
            }
            ystar <- trimr(ystar, max(max_q_F, max_l_F), 0) # (bigT - max) by K_F
            xstar <- trimr(xstar, max(max_q_F, max_l_F), 0) # (bigT - max) by K_F
            TT <- nrow(ystar)
            
            # bring to right dimension
            prior_mean <- matrix(getElement(prior$Psi, "mean"), nrow = q_Fk, ncol = 1) # q_Fk by 1
            inv_prior_var <- diag(1 / getElement(prior$Psi, "var"), nrow = q_Fk) # q_Fk by q_Fk
            sig_inv <- 1/sig2_F[k] # scalar
            
            # posterior mean and variance
            # Dimension of post_var: either q_Fk by q_Fk or K_F by K_F
            # Dimension of post_mean: either 
            post_var <- solve( inv_prior_var + (sig_inv * t(xstar) %*% xstar) ) # 
            post_mean <- post_var %*% (inv_prior_var %*% prior_mean + sig_inv * t(xstar) %*% ystar)
            
            # draw from multivariate normal with posterior mean and variance
            C <- chol(sig2_F[k] * post_var)
            accept <- 0; counter <- 0
            while (accept == 0 & counter < 1e3) {
                psi_F1 <- post_mean + t(C) * rnorm(q_Fk, 1)
                a <- rbind(1, -psi_F1); ceof <- a[nrow(a):1, ]
                root <- polyroot(ceof)
                rootmod <- abs(root)
                if (min(rootmod) > 1.0001) {
                    accept = 1
                } else {
                    counter <- counter + 1
                    accept <- 0
                }
            }
            if (!accept && counter >= 1000) {
                print(counter)
            }
            SSE <- t(ystar - xstar %*% psi_F1) %*% (ystar - xstar %*% psi_F1)
        } else {
            TT <- nrow(Ft)
            ystar <- as.matrix(Ft[, k])
            SSE <- t(ystar) %*% ystar
        }
        
        if (!is.na(var_fix)) {
            sig2_F1 <- var_fix
        } else {
            d = getElement(prior$Sigma, "shape") + SSE
            c = rchisq(n = 1, df = TT + getElement(prior$Sigma, "dof"))
            sig2_F1 <- d/c
        }
        
        xpsi_F1[1:q_Fk] <- psi_F1[1:q_Fk]
        
        if (k == 1) {
            SIG <- sig2_F1
            PSI <- t(xpsi_F1)
        } else {
            SIG <- rbind(SIG, sig2_F1)
            PSI <- rbind(PSI, t(xpsi_F1))
        }
        
        output$SIG <- SIG 
        output$PSI <- PSI
        return(output)
    }
}