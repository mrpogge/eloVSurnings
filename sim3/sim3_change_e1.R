#-------------------------------------------------------------------------------
# Simulation study with changing abilities using the Urnings algorithm.
#-------------------------------------------------------------------------------

#Dependencies
source("algo/elo_wrapper.R")
source("analysis/sim3_settings.R")
library(tictoc)
#-------------------------------------------------------------------------------

#select the first adaptivity-ability distribution condition

Theta = Theta_array[,,1]
mu_p = Mu_p[1]

#container
res_elo = matrix(0, ncol = 4, nrow = length(K_s))
res_bias = res_var =  matrix(0, ncol = 2 + (n_games * n_students), nrow = length(K_s))
colnames(res_elo) = c("adapt_dist", "step_size", "bias", "var")
res_elo[,1] = res_bias[,1] = res_var[,1] = 1
res_elo[,2] = res_bias[,2] = res_var[,2] = K_s

for(us in 1:length(K_s)){
  tic()
  
  print(paste0("Remaining iterations: ", length(K_s)-us))
  
  k_s = K_s[us]
  r = eloChange(n_students = n_students,
                n_items = n_items,
                n_reps = 20,
                n_games = n_games,
                K_s = k_s,
                K_i = K_i,
                Theta = Theta,
                Delta = Delta,
                m_p = mu_p,
                s_p = sigma_p,
                OS = "MAC")
  
  res_elo[us, "bias"] = mean(r$bias_v)
  res_elo[us, "var"] = mean(r$var_v)
  res_bias[us, 3:(n_games * n_students + 2)] = r$bias_v
  res_var[us, 3:(n_games * n_students + 2)] = r$var_v
  toc()
}

#------------------------------------------------------------------------------
#saving the output

saveRDS(res_elo, "output/res3_e1.RDS")
saveRDS(res_bias, "output/res3_bias_e1.RDS")
saveRDS(res_var, "output/res3_var_e1.RDS")
