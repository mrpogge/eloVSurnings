#-------------------------------------------------------------------------------
#                 Baseline calculations for Elo with N(0,1) dist
# NOTE: if the abilties or the urn sizes that are used in the simulations change
# this matrix needs to be updated.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Dependencies
#-------------------------------------------------------------------------------
source("algo/baseline_functions.R")
source("sim1/sim1_settings.R")
source("util.R")

library(tictoc)
library(tidyverse)
#-------------------------------------------------------------------------------
# getting baseline MSE for ELO 
Theta = Theta_matrix[,1]
mu_p = Mu_p[1]
# initialise containers
b_sys_elo = matrix(0, ncol = 2, nrow = length(K_s))
b_elo = matrix(0, ncol = 7, nrow = length(K_s) * length(Theta))
colnames(b_sys_elo) = c("K", "b_mse_sys")
colnames(b_elo) = c("K_s", "ability", "b_mse", "lq", "uq", "meanT", "varT")


# run all K systems
for(nk in 1:length(K_s)){
  tic()
  k = K_s[nk]
  print(paste0("Running K = ", k))
  res = baselineSimpleElo(n_students = n_students,
                          n_items = n_items,
                          n_reps = n_reps_baseline,
                          n_games = n_games_baseline,
                          K_s = k,
                          Theta = Theta,
                          Delta = Delta)
  output = matrix(0,nrow = n_students, ncol = 7)
  output[,1] = k
  output[,2] = Theta
  output[,3] = get_baseline(res$MSE_i, window_size = window_size)
  meanT = get_baseline(res$M_i, window_size = window_size)
  varT = get_baseline(res$Var_i, window_size = window_size)
  output[,4] = meanT - qnorm(0.975) * varT/n_reps #change the value to match the simulation replications
  output[,5] = meanT + qnorm(0.975) * varT/n_reps
  output[,6] = meanT
  output[,7] = varT
  b_elo[((nk-1)*n_students + 1):(nk*n_students), ] = output
  b_sys_elo[nk,] = c(k, mean(output[,3]))
  toc()
}

b_elo = b_elo %>% as.data.frame()
b_sys_elo = b_sys_elo %>% as.data.frame()

#-------------------------------------------------------------------------------
# saving the output (RDS)
saveRDS(b_elo, "output/b1_e1.RDS")
saveRDS(b_sys_elo, "output/b1_sys_e1.RDS")

