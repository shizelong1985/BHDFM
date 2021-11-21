bigZ <- bigZ_ns

paramsF <- list()
paramsF$B <- B
paramsF$q_G <- q_G
paramsF$K_blocks <- K_G
paramsF$q_F <- q_F*matrix(1, K_F, 1)
paramsF$K_facs <- K_F
paramsF$l_F <- l_F

# Initialize factors as PCA estimates plus random noise
g <- vector(mode = "list", length = B)
h <- vector(mode = "list", length = sum(S_b != 0))

for (b in 1:B) {
    g[[b]] <- G_pc[[b]] + noise_var*rnorm(bigT)
    
    if (S_b[b] == 0) { # S_b[b] equals zero if there is no sub-block structure
        bigZ[[b]] <- scale(bigZ[[b]])
        for (k in 1:K_G[b]) {
            g[[b]][,k] <- sign(cor(g[[b]][,k], bigZ[[b]])) * g[[b]][,k]
        }
    } else { # 
        for (s in 1:S_b[b]) {
            bigZ[[b]][[s]] <- scale(bigZ[[b]][[s]])
            h[[b]][[s]] <- H_pc[[b]][[s]] + noise_var*rnorm(bigT)
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
f <- pc_T(allG, K_F)$F.hat + noise_var*rnorm(bigT)

for (k in 1:K_F) {
    f[, k] <- sign(cor(f[, k], allG[, k])) * f[, k]
}

startprior <- list()
startprior$Psi <- list(mean = 0, var = 1e10)
startprior$Sigma <- list(shape = 0, dof = 0)

tmp <- draw_psi_sigma(F_pc, 0.5*matrix(1, nrow = K_F, ncol = 1), matrix(1, nrow = K_F, 1), startprior, paramsF, sig_fix)
