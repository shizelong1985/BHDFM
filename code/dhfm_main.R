
#---------------------------------------------------------------------------------------------------
#
#  Replication code for the paper 'Dynamic Hierarchical Factor Models'
#  R translation created by Anthony M. Thomas, Jr. (amthomasjr@gmail.com)
#  Implements a Gibbs sampler to estimate the following factor model:
#
#       Z_bsnt = Lambda_Hbsn(L)'*H_bst + e_Zbsnt
#       H_bst = Lambda_Gbs(L)'*G_bt + e_Hbst
#       G_bt = Lambda_Fb(L)'*F_t + e_Gbt
#       e_Zbsnt = psi_Zbsn(L)*e_Zbsn_{t-1} + epsilon_Zbsnt
#       e_Hbst = psi_Hbs(L)*e_Hbs_{t-1} + epsilon_Hbst
#       e_Gbt = psi_Gb(L)*e_Gb_{t-1} + epsilon_Gbt
#       F_t = psi_F*F_{t-1} + epsilon_Ft
#       
#---------------------------------------------------------------------------------------------------

n_burn  <- 50000
n_keep  <- 50000
n_skip  <- 50
n_gibbs <- n_burn + n_keep

sig_fix     <- .1
noise_var   <- .2

