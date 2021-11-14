#' ---
#' ---
#' title: Replication code for the paper 'Dynamic Hierarchical Factor Models'
#' author: Anthony M. Thomas Jr.
#' copyright: Emanuel Moench (emanuel.moench@ny.frb.org), Serena Ng (serena.ng@columbia.edu), Simon Potter (simon.potter@ny.frb.org) 
#' ---


library(bannerCommenter) # to create comment banners
library(R.matlab)
library(matlab)
library(mSTEM) # for 'conv' function
library(signal) # for 'filter' function
library(fBasics) # for 'vec' function

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##                                                                       %%
##  Replication code for the paper 'Dynamic Hierarchical Factor Models'  %%
##                                                                       %%
##                     %%
##                                                                       %%
##  R translation created by Anthony Thomas (amthomasjr@gmail.com)       %%
##                                                                       %%
##  Implements a Gibbs sampler to estimate the following factor model:   %%
##                                                                       %%
##       Z_bsnt = Lambda_Hbsn(L)'*H_bst + e_Zbsnt                        %%
##       H_bst = Lambda_Gbs(L)'*G_bt + e_Hbst                            %%
##       G_bt = Lambda_Fb(L)'*F_t + e_Gbt                                %%
##       e_Zbsnt = psi_Zbsn(L)*e_Zbsn_{t-1} + epsilon_Zbsnt              %%
##       e_Hbst = psi_Hbs(L)*e_Hbs_{t-1} + epsilon_Hbst                  %%
##       e_Gbt = psi_Gb(L)*e_Gb_{t-1} + epsilon_Gbt                      %%
##       F_t = psi_F*F_{t-1} + epsilon_Ft                                %%
##                                                                       %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Intro ----
n_burn <- 50000; n_keep <- 50000; n_skip <- 50; n_gibbs <- n_burn + n_keep;

## Load data ----

matlab_data <- readMat('~/OneDrive\ -\ University\ of\ Texas\ at\ Arlington/Research\ Project/References/Dynamic\ Hierarchical\ Factor\ Models/dataverse_files/dhfm_data.mat')
data <- matlab_data$bigZ.ns

## Define model parameters ----

B <- 5 # number of blocks
K_F <- 1 # number of common factors
K_G <- c(1, 1, 1, 1, 1) # K_G[b] is number of block-level factors in block b
K_H <- vector(mode = "list", length = B) # K_H[[b]] the number of subblock-level factors in block b
K_H[[1]] <- c(1, 1, 2); K_H[[2]] <- c(2, 1); K_H[[3]] <- c(1, 1)
K_H[[4]] <- 0; K_H[[5]] <- 0

j # makes the hierarchy more logical due to import oddness from MATLAB
for (b in 1:B) {
        data[[b]] <- data[[b]][[1]]
        if (is.list(data[[b]])) {
                for (s in 1:length(K_H[[b]])) {
                        data[[b]][[s]] <- data[[b]][[s]][[1]]
                }
        }
}
bigZ_ns <- data
bigZ <- bigZ_ns
bigT <- 227

Bsub <- rep(NA, B) # Bsub[b] is the number of subblocks in block b
for (b in 1:B) {
        Bsub[b] <- length(K_H[[b]])
        if (length(K_H[[b]]) == 1 && unlist(K_H[[b]]) == 0) {
                Bsub[b] <- 0
        }
}

q_F <- 1 # lag order of transition equation
q_G <- rep(1, B)
q_H <- rep(1, B)
q_Z <- rep(1, B)

l_F <- rep(0, B)
l_G <- rep(0, B)
l_H <- rep(0, B)

Nsub <- matrix(0, nrow = sum(Bsub != 0), ncol = sum(Bsub != 0)) # Nsub[b, s] is the size of subblock s of block b
for (b in 1:sum(Bsub != 0)) {
        if (Bsub[b] != 0) {
                for (s in 1:Bsub[b]) {
                        Nsub[b ,s] <- ncol(data[[b]][[s]])
                }
        }
}

sig_fix = .1;
noise_var = .2;

## Functions ----
pc_T <- function (y, nfac) {
        bigt <- nrow(y); bign <- ncol(y)
        yy <- y %*% t(y)
        svd.yy <- svd(yy)
        Fhat0 <- svd.yy$u; eigval <- diag(svd.yy$d); Fhat1 <- svd.yy$v
        fhat <- Fhat0[, 1:nfac] * sqrt(bigt)
        lambda <- t(y) %*% fhat / bigt
        ehat <- y - fhat %*% t(lambda)
        ss <- diag(eigval)
        returnList <- list(lambda = lambda, fhat = fhat, ehat = ehat, ss = ss)
}
trimr <- function (x, a, b) {
        nt <- nrow(x)
        nc <- ncol(x)
        if (a >= 0) {
                xx <- as.matrix(x[(a + 1):nt, ])
        }
        if (b >= 0) {
                if (a > 0) {
                        x <- xx
                }
                nt <- nrow(x)
                nc <- ncol(x)
                xx <- as.matrix(x[1:(nt - b), 1:nc])
        }
        return(xx)
}
lagn <- function (x, n) {
        nt <- nrow(x)
        nc <- ncol(x)
        if (n > 0) {
                x1 <- trimr(x, 0, n)
                y <- rbind(matrix(0, nrow = n, ncol = nc), x1)
        }
        if (n < 0) {
                x1 <- trimr(x, abs(n), 0)
                y <- rbind(x1, matrix(0, nrow = abs(n), ncol = nc))
        }
        if (n == 0) {
                y <- x
        }
        return(y)
}
draw_psi_sigma <- function (Ft, psi_F, sig2_F, prior.psi.mean, prior.psi.var, prior.Sigma.shape, prior.Sigma.dof, l_F, q_F, var_fix) {
        K_F <- ncol(Ft)
        max_q_F <- max(q_F)
        max_l_F <- max(l_F)
        bigT <- nrow(Ft)
        PSI <- psi_F; SIG <- sig2_F
        
        for (k in 1:K_F) {
                q_Fk <- q_F[k]
                xpsi_F1 = matrix(0, nrow = max_q_F, 1)
                if (q_Fk > 0) {
                        ystar <- as.matrix(Ft[, k])
                        xstar <- matrix(, nrow = bigT, ncol = 0)
                        for (i in 1:q_Fk) {
                                xstar <- cbind(xstar, lagn(f, i))
                        }
                        ystar <- trimr(ystar, max(max_q_F, max_l_F), 0)
                        xstar <- trimr(xstar, max(max_q_F, max_l_F), 0)
                        TT <- nrow(ystar)
                        
                        prior_mean <- prior.psi.mean * matrix(1, nrow = q_Fk, ncol = 1)
                        inv_prior_var <- as.numeric(solve(prior.psi.var)) * diag(1, nrow = q_Fk)
                        sig_inv <- 1/sig2_F[k]
                        
                        post_var <- solve(inv_prior_var + sig_inv * t(xstar) %*% xstar)
                        post_mean <- post_var %*% (inv_prior_var %*% prior_mean + sig_inv %*% t(xstar) %*% ystar)
                        
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
                        # if (!accept & counter >= 1000) {
                        #      print(counter)
                        # }
                        SSE <- t(ystar - xstar %*% psi_F1) %*% (ystar - xstar %*% psi_F1)
                } else {
                        TT <- nrow(Ft)
                        ystar <- as.matrix(Ft[, k])
                        SSE <- t(ystar) %*% ystar
                }
                if (!is.null(var_fix)) {
                        sig2_F1 <- var_fix
                } else {
                        d <- prior.Sigma.shape + SSE
                        c <- rchisq(1, df = TT + prior.Sigma.dof)
                        sig2_F1 <- d/c
                }
                xpsi_F1[1:q_Fk] <- psi_F1[1:q_Fk]
                if (k == 1) {
                        SIG <- sig2_F1
                        PSI <- xpsi_F1
                } else {
                        SIG <- rbind(SIG, sig2_F1)
                        PSI <- rbind(PSI, t(xpsi_F1))
                }
        }
        output <- list(SIG = SIG, PSI = PSI)
        return(output)
}
compute_resids <- function (obs, fac, Lambda, K_F, l_F) {
        max_l_F <- max(l_F)
        Fmat <- fac
        fit <- matrix(0, nrow = nrow(obs), ncol = ncol(obs))
        fit <- fit + Fmat %*% t(Lambda)
        e <- obs - fit
        return(e)
}
switch_rows <- function (X, id) {
        id1 <- id; id2 <- flipud(id)
        temp <- X
        X[id1, ] <- temp[id2, ]
        X[id2, ] <- temp[id1, ]
        return(X)
}
switch_cols <- function (X, id) {
        id1 <- id; id2 <- t(flipud(t(id)))
        temp <- X
        X[, id1] <- temp[, id2]
        X[, id2] <- temp[, id1]
        return(X)
}

# Initialize Sampler ----

## Preliminary PC analysis ----

H_pc <- vector(mode = "list", length = 3)
H_pc[[1]] <- list(pc_T(scale(data[[1]][[1]]), K_H[[1]][1])$fhat, 
                  pc_T(scale(data[[1]][[2]]), K_H[[1]][2])$fhat, 
                  pc_T(scale(data[[1]][[3]]), K_H[[1]][3])$fhat)
H_pc[[2]] <- list(pc_T(scale(data[[2]][[1]]), K_H[[2]][1])$fhat, 
                  pc_T(scale(data[[2]][[2]]), K_H[[2]][2])$fhat)
H_pc[[3]] <- list(pc_T(scale(data[[3]][[1]]), K_H[[3]][1])$fhat, 
                  pc_T(scale(data[[3]][[2]]), K_H[[3]][2])$fhat)

G_pc <- vector(mode = "list", length = B)
G_pc[[1]] <- pc_T(scale(matrix(unlist(H_pc[[1]]), nrow = bigT)), K_G[1])$fhat
G_pc[[2]] <- pc_T(scale(matrix(unlist(H_pc[[2]]), nrow = bigT)), K_G[2])$fhat
G_pc[[3]] <- pc_T(scale(matrix(unlist(H_pc[[3]]), nrow = bigT)), K_G[3])$fhat
G_pc[[4]] <- pc_T(scale(data[[4]]), K_G[4])$fhat
G_pc[[5]] <- pc_T(scale(data[[5]]), K_G[5])$fhat

F_pc <- pc_T(matrix(unlist(G_pc), nrow = bigT), K_F)$fhat

GX_pc <- vector(mode = "list", length = B)
for (b in 1:B) {
        GX_pc[[b]] <- pc_T(scale(matrix(unlist(data[[b]]), nrow = bigT)), K_G[b])$fhat
}
FX_pc <- pc_T(scale(matrix(unlist(data), nrow = bigT)), K_F)$fhat
FH_pc <- pc_T(scale(matrix(unlist(H_pc), nrow = bigT)), 8)$fhat

## Initialize factors as PC estimates plus random noise ----

g <- vector(mode = "list", length = B)
h <- vector(mode = "list", length = sum(Bsub != 0))

for (b in 1:B) {
        g[[b]] <- matrix(G_pc[[b]][1:K_G[b]] + noise_var*rnorm(bigT), nrow = bigT)
        if (Bsub[b] == 0) {
                bigZ[[b]] <- scale(bigZ[[b]])
                for (k in 1:K_G[b]) {
                        g[[b]][, k] <- sign(cor(g[[b]][, k], bigZ[[b]][, k])) * g[[b]][, k]
                }
        } else {
                for (s in 1:Bsub[b]) {
                        bigZ[[b]][[s]] <- scale(bigZ[[b]][[s]])
                        h[[b]][[s]] <- matrix(H_pc[[b]][[s]] + noise_var*rnorm(bigT), nrow = bigT)
                        for (k in 1:K_H[[b]][s]) {
                                h[[b]][[s]][, k] <- sign(cor(h[[b]][[s]][, k], bigZ[[b]][[s]][, k])) * h[[b]][[s]][, k]
                        }
                }
                allH <- matrix(unlist(h[[b]]), nrow = bigT)
                for (k in 1:K_G[b]) {
                        g[[b]][, k] <- sign(cor(g[[b]][, k], allH[, k])) * g[[b]][, k]
                }
        }
}

allG <- matrix(unlist(g), nrow = bigT)
f <- matrix(pc_T(allG, K_F)$fhat, bigT) + noise_var*rnorm(bigT)

for (k in 1:K_F) {
        f[, k] <- sign(cor(f[, k], allG[, k])) * f[, k]
}

## Initialize parameters based on starting values for the factors ----
# First regress starting values of aggregate factors on their own lags to obtain starting values for the psi_F's and sig2_F's

startprior <- list(Psi = list(mean = 0, var = 1e10), Sigma = list(shape = 0, dof = 0))
args <- list(Ft = f, psi_F = 0.5*matrix(1, nrow = K_F, ncol = 1), 
             sig2_F = matrix(1, nrow = K_F, ncol = 1), 
             prior.psi.mean = startprior$Psi$mean, prior.psi.var = startprior$Psi$var, 
             prior.Sigma.shape = startprior$Sigma$shape, prior.Sigma.dof = startprior$Sigma$dof, 
             l_F = l_F, q_F = q_F, var_fix = sig_fix)
psi_F <- do.call(draw_psi_sigma, args)$PSI; sig2_F <- do.call(draw_psi_sigma, args)$SIG

F_lags <- f
lambda_F <- t( solve(t(F_lags) %*% F_lags, t(F_lags) %*% allG) )
Lambda_F <- as.matrix(lambda_F[, 1])

signlam <- NULL
for (i in 1:K_F) {
     signlam[i] <- sign(Lambda_F[i])
     f[, i] <- f[,i] * signlam
     Lambda_F[, i] <- Lambda_F[, i] * signlam[i]
}
Lambda_F[1, 1] <- 1

eG <- vector(mode = "list", length = B)
eZ <- vector(mode = "list", length = B)
eH <- vector(mode = "list", length = B)
psi_G <-vector(mode = "list", length = B)
psi_Z <- vector(mode = "list", length = B)
psi_H <- vector(mode = "list", length = sum(Bsub != 0))
sig2_G <- vector(mode = "list", length = B)
sig2_Z <- vector(mode = "list", length = B)
sig2_H <- vector(mode = "list", length = B)
Lambda_G <- vector(mode = "list", length = B)
Lambda_Gsub <- vector(mode = "list", length = B)

Lambda_H <- vector(mode = "list", length = sum(Bsub != 0))
for (i in 1:sum(Bsub != 0)) {
        Lambda_H[[i]] <- vector(mode = "list", length = Bsub[i])
        for (j in 1:Bsub[i]) {
                Lambda_H[[i]][[j]] <- matrix(NA, nrow = ncol(data[[i]][[j]]), ncol = K_H[[i]][[j]])
        }
}




for (b in 1:B) {
     lambda_F_mat <- Lambda_F
     lambda_Fb <- as.matrix(lambda_F_mat[b, ])
     eG[[b]] <- compute_resids(g[[b]], f, lambda_Fb, K_F = K_G[b], l_F = l_G[b])
     args.G <- list(Ft = eG[[b]], psi_F = 0.5*matrix(1, nrow = K_G[b], ncol = 1), 
                  sig2_F = matrix(1, nrow = K_G[b], ncol = 1), 
                  prior.psi.mean = startprior$Psi$mean, prior.psi.var = startprior$Psi$var, 
                  prior.Sigma.shape = startprior$Sigma$shape, prior.Sigma.dof = startprior$Sigma$dof, 
                  l_F = l_G[b], q_F = q_G[b]*rep(1, K_G[b]), var_fix = sig_fix)
     psi_G[[b]] <- do.call(draw_psi_sigma, args.G)$PSI
     sig2_G[[b]] <- do.call(draw_psi_sigma, args.G)$SIG
     
     if (Bsub[b] == 0) {
          yy <- bigZ[[b]]
          
          # regress all Z_bn on G_b in order to get starting values for Lambda_Gb
          gb_lags <- g[[b]]
          
          # regress individual variables on subblock factor starting values
          lambda_Gb <- t(solve( t(gb_lags) %*% gb_lags, t(gb_lags) %*% yy ) )
          Lambda_G[[b]] <- as.matrix(lambda_Gb[, 1])

          # Find Zbi whose loading on Gb is closest to 1 in absolute value
          mindist <- apply(abs(Lambda_G[[b]]) - matrix(1, nrow = nrow(Lambda_G[[b]]), ncol = ncol(Lambda_G[[b]])), 2, function(x) {sort(x, decreasing = T)})
          neworder <- apply(abs(Lambda_G[[b]]) - matrix(1, nrow = nrow(Lambda_G[[b]]), ncol = ncol(Lambda_G[[b]])), 2, function(x) {order(x, decreasing = T)})
          Lambda_G[[b]] <- switch_rows(Lambda_G[[b]], neworder[1, ])
          bigZ[[b]] <- switch_cols(bigZ[[b]], neworder[1, ])
          bigZ_ns[[b]] <- switch_cols(bigZ_ns[[b]], neworder[1, ])
          signlam <- sign(Lambda_G[[b]][1, 1])
          g[[b]][, 1] <- g[[b]] * signlam
          Lambda_G[[b]][, 1] <- Lambda_G[[b]][, 1] * signlam
          Lambda_G[[b]][1, 1] <- 1
          
          eZ[[b]] = compute_resids(yy, g[[b]], Lambda = Lambda_G[[b]], K_F = K_G[b], l_F = l_G[b]*matrix(1, nrow = nrow(yy), ncol = 1))
          args.Z <- list(Ft = eZ[[b]], psi_F = 0.5*matrix(1, nrow = K_G[b], ncol = 1), sig2_F = matrix(1, nrow = nrow(yy), ncol = 1), 
                       prior.psi.mean = startprior$Psi$mean, prior.psi.var = startprior$Psi$var, 
                       prior.Sigma.shape = startprior$Sigma$shape, prior.Sigma.dof = startprior$Sigma$dof, 
                       l_F = l_F[b] * matrix(1, nrow = nrow(yy), ncol = 1), q_F = q_Z[b]*matrix(1, nrow = nrow(yy), ncol = 1), var_fix = sig_fix)
          psi_Z[[b]] <- do.call(draw_psi_sigma, args.Z)$PSI; sig2_Z[[b]] <- do.call(draw_psi_sigma, args.Z)$SIG
          
     } else { # if Bsub[b] != 0)
          # Lambda_G <- matrix(, nrow = bigT, ncol = 0)
          # g[[b]] <- g[[b]]
          gb_lags <- g[[b]]
          Hb <- matrix(unlist(h[[b]]), nrow = bigT)

          # regress individual variables on subblock factor starting values
          lambda_Gb <- t( solve(t(gb_lags) %*% gb_lags, t(gb_lags) %*% Hb) )
          Lambda_Gsub[[b]] <- as.matrix(lambda_Gb[, K_G[b]])
          signlam <- sign(Lambda_Gsub[[b]][1, 1])
          g[[b]][, 1] <- g[[b]] * signlam
          Lambda_Gsub[[b]][, 1] <- as.matrix(Lambda_Gsub[[b]][, 1]) * signlam
          Lambda_Gsub[[b]][1, 1] <- 1

          for (s in 1:Bsub[b]) {
               yy <- bigZ[[b]][[s]]
               if (s == 1) {
                    Hbs_indx <- 1:K_H[[b]][s]
               } else {
                    Hbs_indx <- (sum(K_H[[b]][1:(s-1)]) + 1):sum(K_H[[b]][1:s])
               }
               
               lambda_Gsub_mat <- Lambda_Gsub[[b]]
               lambda_Gbs <- lambda_Gsub_mat[Hbs_indx, ]
               eH[[b]][[s]] <- compute_resids(h[[b]][[s]], g[[b]], lambda_Gbs, K_F = K_H[[b]][[s]], l_F = l_H[b]*ones(nrow(h[[b]][[s]]), 1))
               args.H <- list(Ft = eH[[b]][[s]], psi_F = 0.5*matrix(1, nrow = K_H[[b]][[s]], ncol = 1), 
                              sig2_F = matrix(1, nrow = K_H[[b]][[s]], ncol = 1), 
                              prior.psi.mean = startprior$Psi$mean, prior.psi.var = startprior$Psi$var, 
                              prior.Sigma.shape = startprior$Sigma$shape, prior.Sigma.dof = startprior$Sigma$dof, 
                              l_F = l_H[b]*ones(nrow(h[[b]][[s]]), 1), q_F = q_Z[b]*ones(nrow(yy), 1), var_fix = sig_fix)
               psi_H[[b]][[s]] <- do.call(draw_psi_sigma, args.H)$PSI; sig2_H[[b]][[s]] <- do.call(draw_psi_sigma, args.H)$SIG
               
               hbs_lags <- h[[b]][[s]]
               
               lambda_Hbs <- t( solve( t(hbs_lags) %*% hbs_lags, t(hbs_lags) %*% yy ) )
               Lambda_H[[b]][[s]] <- as.matrix(lambda_Hbs[, 1:K_H[[b]][s]])
               
               mindist <- apply(abs(Lambda_H[[b]][[s]]) - ones(nrow(Lambda_H[[b]][[s]]), ncol(Lambda_H[[b]][[s]])), 2, function(x) {sort(x, decreasing = T)})
               neworder <- apply(abs(Lambda_H[[b]][[s]]) - ones(nrow(Lambda_H[[b]][[s]]), ncol(Lambda_H[[b]][[s]])), 2, function(x) {order(x, decreasing = T)})
               
               Lambda_H[[b]][[s]] <- switch_rows(Lambda_H[[b]][[s]], neworder[1, ])
               
               bigZ[[b]][[s]] <- switch_cols(bigZ[[b]][[s]], neworder[1, ])
               bigZ_ns[[b]][[s]] <- switch_cols(bigZ_ns[[b]][[s]], neworder[1, ])
               
               for (i in 1:K_H[[b]][[s]]) {
                       signlam[i] <- sign(Lambda_H[[b]][[s]][i, i])
                       h[[b]][[s]][ , i] <- h[[b]][[s]][, i] * signlam[i]
                       Lambda_H[[b]][[s]][ , i] <- Lambda_H[[b]][[s]][, i] * signlam[i]
                       }
               
               diag(Lambda_H[[b]][[s]]) <- rep(1, ncol(Lambda_H[[b]][[s]]))
               Lambda_H[[b]][[s]][upper.tri(Lambda_H[[b]][[s]])] <- 0
               
               eZ[[b]][[s]] <- compute_resids(yy, h[[b]][[s]], Lambda_H[[b]][[s]], K_F = K_H[[b]][[s]], l_F = l_H[b]*ones(nrow(h[[b]][[s]]), 1))
               args.Z <- list(Ft = eZ[[b]][[s]], psi_F = 0.5*ones(Nsub[b, s], 1), sig2_F = ones(Nsub[b, s], 1), 
                              prior.psi.mean = startprior$Psi$mean, prior.psi.var = startprior$Psi$var, 
                              prior.Sigma.shape = startprior$Sigma$shape, prior.Sigma.dof = startprior$Sigma$dof, 
                              l_F = l_H[b]*ones(nrow(yy), 1), q_F = q_Z[b]*ones(nrow(yy), 1), var_fix = sig_fix)
               psi_Z[[b]][[s]] <- do.call(draw_psi_sigma, args.Z)$PSI; sig2_Z[[b]][[s]] <- do.call(draw_psi_sigma, args.Z)$SIG
               } # end for (s in 1:Bsub[b])
        } # end else if (Bsub[b] == 0)
} # end for (b in 1:B)










# Run Sampler ----

jj <- 1

## Priors ----
prior.F <- list(Lambda.mean = 0, Lambda.var = 10, Psi.mean = 0.5, Psi.var = 10, Sigma.shape = 0.01, Sigma.dof = 4)
prior.G <- list(Lambda.mean = 0, Lambda.var = 10, Psi.mean = 0.5, Psi.var = 10, Sigma.shape = 0.01, Sigma.dof = 4)
prior.H <- list(Lambda.mean = 0, Lambda.var = 10, Psi.mean = 0.5, Psi.var = 10, Sigma.shape = 0.01, Sigma.dof = 4)
prior.Z <- list(Psi.mean = 0.5, Psi.var = 10, Sigma.shape = 0.01, Sigma.dof = 4)
priors <- list(prior.F, prior.G, prior.H, prior.Z)

cal_alpha <- function(FF, Psi, Lambda) {
        
        K_F <- size(FF, 2)
        l_F <- length(Lambda_F) - 1
        q_G <- size(Psi, 2)
        Lambda_lags <- NULL
        for (l in 1:(l_F + 1)) {
                Lambda_lags <- cbind(Lambda_lags, Lambda[[l]])
        }
        l_star <- l_F + q_G
        
        B <- size(Lambda[[1]], 1)
        Lambda_tilde <- NULL
        alpha <- NULL
        for (b in 1:B) {
                lambdab_tilde <- NULL
                for (j in 1:K_F) {
                        indx <- K_F:(K_F*(l_F + 1))
                        lambdabj_tilde <- conv(as.vector(cbind(1, -Psi[b, ])), as.vector(Lambda_lags[b, indx]))
                        lambdab_tilde <- cbind(lambdab_tilde, t(lambdabj_tilde))
                }
                lambdab_tilde <- reshape(t(lambdab_tilde), 1, K_F*(l_star + 1))
                Lambda_tilde <- rbind(Lambda_tilde, lambdab_tilde)
                Flags <- NULL
                for (l in 0:l_star) {
                        Flags <- cbind(Flags, lagn(FF, l))
                }
                alpha[[b]] <- Flags %*% t(lambdab_tilde)
        }
        output <- list(alpha, Lambda_tilde)
        return(output)
}

sample_facs <- function(g, alpha, lambda_tilde, psi_F, sig2_F, psi_G, sig2_G, params) {
        
        B <- params[["B"]]
        q_G <- params[["q_G"]]
        max_q_G <- max(q_G)
        K_G <- params[["K_G"]]
        K_F <- params[["K_facs"]]
        q_F <- params[["q_F"]]
        max_q_F <- max(q_F)
        
        Gmat <- NULL
        sig2_G_vec <- NULL
        psi_G_mat <- NULL
        
        for (b in 1:B) {
                Gmat <- cbind(Gmat, g[[b]])
                sig2_G_vec <- rbind(sig2_G_vec, sig2_G[[b]])
                
                if (max_q_G > 0) {
                        psi_G_mat <- rbind(psi_G_mat, cbind(psi_G[[b]], matrix(0, nrow = K_G[b], ncol = max_q_G - q_G[b])))
                } else {
                        psi_G_mat <- rbind(psi_G_mat, matrix(0, nrow = K_G[b], ncol = 1))
                }
        }
        TT <- size(Gmat, 1); N <- size(Gmat, 2)
        ystar <- matrix(0, nrow = TT, ncol = N)
        for (i in 1:N) {
                require(signal)
                ystar[ , i] <- filter(as.vector(cbind(1, -Psi[b, ])), 1, Gmat[ , i])
        }
        ystar <- t(ystar)
        
        # transition matrix of the filtered data is r = p + 1 + q + 1
        ndim <- size(lambda_tilde, 2)
        psi_F_mat <- zeros(K_F, max_q_F * K_F)
        for (r in 1:K_F) {
                indx <- seq(r, K_F*max_q_F, K_F)
                psi_F_mat[r, indx] <- psi_F[r, ]
        }
        
        A <- zeros(ndim)
        A[1:size(psi_F_mat, 1), 1:(ndim - K_F)] <- eye(ndim - K_F)
        H <- lambda_tilde
        QQ <- zeros(ndim, ndim)
        QQ[1:K_F, 1:K_F] <- diag(sig2_F)
        R <- diag(sig2_G_vec)
        Alpha <- cbind(alpha, zeros(nrows(alpha), ndim - ncol(alpha)))
        
        # Initialize Kalman Filter
        Gtt <- zeros(ndim, 1)
        P00 <- solve(eye(ndim^2) - kronecker(A, A)) %*% vec(QQ)
        Ptt <- reshape(P00, ndim, ndim)
        
        # Kalman filter recursion
        Fmat <- zeros(ndim, TT)
        Pmat <- vector(mode = "list", length = TT)
        
        for (t in 1:TT) {
                Gtt1    <- t(Alpha[t, ]) + A %*% Gtt1           # x(t|t-1)
                Ptt1    <- A %*% Ptt %*% t(A) + QQ              # P(t|t-1)
                ett1    <- ystar[ , t] - H %*% Gtt1             # eta(t|t-1)=y(t)- y(t|t-1)= y(t)-H*x(t|t-1)
                v_ett1  <- H %*% Ptt11 %*% t(H) + R             # var(eta(t|t-1))
                k_gn    <- Ptt1 %*% t(H) %*% solve(v_ett1)      # K=P(t|t-1)H'/ f(t|t-1)
                Gtt     <- Gtt1 + k_gn %*% ett1                 # x(t|t)=x(t|t-1)+ K eta(t|t-1)
                Ptt     <- (eye(ndim) - k_gn %*% H) %*% Ptt1    # P(t|t)= P(t|t-1)-K*H*P(t|t-1)
                
                Fmat[ , t] <- Gtt
                Pmat[[t]] <- Ptt
        }
        
        # Carter-Kohn backward sampling algorithm
        FT_mat <- zeros(TT, K_F)
        G <- t(chol(Pmat[[TT]]))
        FT <- Fmat[ , TT] + G %*% matrix(rnorm(ndim), nrow = ndim, ncol = 1)
        FT_mat[TT, 1:K_F] <- FT[1:K_F, 1]
        jj <- 1:K_F
        
        for (t in (TT-1):1) {
                etT <- FT[jj] - t(Alpha[t + 1, jj]) - A[jj, ] %*% Fmat[ , t]
                v_etT <- A[jj, ] %*% Pmat[[t]] %*%t(A[jj, ]) + QQ[jj, jj]
                k_gn0 <- Pmat[[t]] %*% t(A[jj, ]) %*% solve(v_etT)
                Fmat[ , t] <- Fmat[ , t] + k_gn0 %*% etT
                Pmat[[t]] <- (eye(ndim) - k_gn0 %*% A[jj, ]) %*% Pmat[[t]]
                G <- chol(Pmat[[t]])
                FT <- Fmat[ , t] + G %*% rnorm(ndim)
                FT_mat[t, 1:K_F] <- FT[1:KF, 1]
        }
        done <- 1
        return(FF_mat)
}

draw_lambda <- function(G, FF, psi_G, sig2_G, prior, params) {
      
       B <- params[["B"]]
       q_G <- params[["q_G"]]
       max_q_G <- max(q_G)
       K_G <- params[["K_blocks"]]
       l_F <- params[["l_F"]]
       max_l_F <- max(l_F)
       Lambda_F <- NULL
       count <- 0
       LAMBDA <- vector(mode = "list", length = max_l_F + 1)
       
       for (b in 1:B) {
               Gb <- G[[b]]
               TT <- size(Gb, 1)
               k_G <- size(Gb, 2)
               l_Fb <- l_F[b]
               
               for (i in 1:k_G) {
                       count <- count + 1
                       
                       F_filt <- filter(cbind(1, -psi_G[[b]][i, ]), 1, FF)
                       ystar <- filter(cbind(1, -psi_G[[b]][i, ]), 1, Gb[ , i])
                       Fstar <- F_filt
                       for (j in 1:l_Fb) {
                               Fstar <- cbind(Fstar, lagn(F_filt, j))
                       }
                       Fstar <- trimr(Fstar, max(max_q_G, l_Fb), 0)
                       ystar <- trimr(ystar, max(max_q_G, l_Fb), 0)
                       
                       # bring to right dimension
                       prior_mean <- prior[["mean"]] %*% ones(K_F*(l_Fb + 1), 1)
                       prior_var <- prior[["var"]] %*% eye(K_F*(l_Fb + 1))
                       sig_inv <- 1/sig2_G[[b]][i]
                       
                       # posterior mean and variance
                       post_var <- solve( solve(prior_var) + sig_inv %*% t(Fstar) %*% Fstar )
                       post_mean <- post_var %*% ( solve(prior_var) %*% prior_mean + sig_inv %*% t(Fstar) %*% ystar )
                       
                       # draw from multivariate normal with posterior mean and variance
                       C <- chol(post_var)
                       lambda_Fi <- post_mean + t(C) %*% matrix(rnorm(K_F * (l_Fb + 1)), ncol = 1)
                       Lambda_F <- rbind(Lambda_F, t(lambda_Fi))
                       
               } # end for (i in 1:k_G)
       } # end for (b in 1:B)
       
       # Impose lower triangularity with ones on the diagonal for identification
       
       for (l in 1:(max_l_F + 1)) {
              LAMBDA[[l]] <- Lambda_F[ , ((l-1)*K_F + 1):(l*K_F)]
              for (j in 1:K_F) {
                if (l == 1) {
                  LAMBDA[[l]][j, j] <- 1
                }
                LAMBDA[[l]][j, (j+1):K_F] <- 0
              }
       }
}

Lambda_F <- list(Lambda_F)
for (gibbs in 1:n_gibbs) {
        # Update global factors and corresponding parameters using last iteration's draws of params and block factors
        alpha_F <- cal_alpha(f, matrix(unlist(psi_G), ncol = 1), Lambda_F)$alpha
        lambda_F_tilde <- cal_alpha(f, unlist(psi_G), Lambda_F)$Lambda_tilde
        
        paramsF <- list(B <- B, q_G <- q_G, q_F <- q_F, K_G <- K_G, K_F <- K_facs)
        f <- do.call(sample_facs, c(g, zeros(size(f)), lambda_F_tilde, psi_F, sig2_F, psi_G, sig2_G, paramsF))
} # end for (gibbs in 1:n_gibbs)


























