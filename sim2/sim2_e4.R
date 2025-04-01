#-------------------------------------------------------------------------------
#                       Simulation 2: adaptivity / dist case 1
#-------------------------------------------------------------------------------
library(tictoc)
# load settings and functions
source("analysis/sim2_settings.R")
source("algo/elo_wrapper.R")

Theta = Theta_matrix[,4]
mu_p = Mu_p[4]
#load baselines
b_mse_sys_elo = readRDS("output/b2_sys_e4.RDS")
b_mse_elo = readRDS("output/b2_e4.RDS")

#initialize containers
res_sys_elo = matrix(0, ncol = 3, nrow = length(K_s))
colnames(res_sys_elo) = c("K", "b_mse_sys", "HT")
res_sys_elo[,1:2] = as.matrix(b_mse_sys_elo)

res_elo = matrix(0, ncol = 4, nrow = length(K_s) * length(Theta))
colnames(res_elo) = c("K", "ability", "b_mse", "HT")
res_elo[,1:3] = as.matrix(b_mse_elo[,1:3])

# run all K systems
for(nk in 1:length(K_s)){
  tic()
  #setting k
  k = K_s[nk]
  #getting baseline for k
  lq_baseline = b_mse_elo[b_mse_elo$K_s == k, 4]
  uq_baseline = b_mse_elo[b_mse_elo$K_s == k, 5]
  
  #status print
  print(paste0("Running K = ", k))
  
  #running the system
  res = eloDoubleHT(n_students = n_students,
                    n_items = n_items,
                    n_reps = n_reps,
                    n_games = n_games,
                    K_s = k,
                    K_i = K_i,
                    Theta = Theta,
                    Delta = Delta,
                    lq_baseline = lq_baseline,
                    uq_baseline = uq_baseline,
                    m_p = 0,
                    s_p = 1,
                    OS = "MAC")
  
  #saving the results to the output files 
  res_elo[b_mse_elo$K == k, "HT"] = res$HT_sys
  res_sys_elo[nk,3] = mean(res$HT_sys)
  toc()
}

#saving the output
saveRDS(res_sys_elo, "output/res2_sys_e4.RDS")
saveRDS(res_elo, "output/res2_e4.RDS")
