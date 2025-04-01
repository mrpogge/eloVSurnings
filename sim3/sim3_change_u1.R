#-------------------------------------------------------------------------------
# Simulation study with changing abilities using the Urnings algorithm.
#-------------------------------------------------------------------------------

#Dependencies
source("algo/urnings_wrapper.R")
source("analysis/sim3_settings.R")
library(tictoc)
#-------------------------------------------------------------------------------

#select the first adaptivity-ability distribution condition

Theta = Theta_array[,,1]
mu_p = Mu_p[1]

#container
res_urnings = matrix(0, ncol = 4, nrow = length(Student_urn_size))
res_bias = res_var =  matrix(0, ncol = 2 + (n_games * n_students), nrow = length(Student_urn_size))
colnames(res_urnings) = c("adapt_dist", "step_size", "bias", "var")
res_urnings[,1] = res_bias[,1] = res_var[,1] = 1
res_urnings[,2] = res_bias[,2] = res_var[,2] = Student_urn_size

for(us in 1:length(Student_urn_size)){
  tic()
  
  print(paste0("Remaining iterations: ", length(Student_urn_size)-us))
  
  student_urn_size = Student_urn_size[us]
  r = urningsChange(n_students = n_students,
                    n_items = n_items,
                    n_reps = 20,
                    n_games = n_games,
                    student_urn_size = student_urn_size,
                    item_urn_size = item_urn_size,
                    Theta = Theta,
                    Delta = Delta,
                    adaptive = 1,
                    m_p = mu_p,
                    s_p = sigma_p,
                    paired = 1,
                    OS = "MAC")
  
  res_urnings[us, "bias"] = mean(r$bias_v)
  res_urnings[us, "var"] = mean(r$var_v)
  res_bias[us, 3:(n_games * n_students + 2)] = r$bias_v
  res_var[us, 3:(n_games * n_students + 2)] = r$var_v
  toc()
}

#------------------------------------------------------------------------------
#saving the output

saveRDS(res_urnings, "output/res3_u1.RDS")
saveRDS(res_bias, "output/res3_bias_u1.RDS")
saveRDS(res_var, "output/res3_var_u1.RDS")
