#-------------------------------------------------------------------------------
#                       Simulation 2: dist/algo case 1
#-------------------------------------------------------------------------------
library(tictoc)
# load settings and functions
source("sim2/sim2_settings.R")
source("algo/urnings_wrapper.R")

Theta = Theta_matrix[,3]
m_p = Mu_p[3]

#load baselines
b_mse_urnings = as.data.frame(readRDS("output/b_u3.RDS"))
b_mse_sys_urnings = as.data.frame(readRDS("output/b_sys_u3.RDS"))

#initialize containers
res_sys_urnings = matrix(0, ncol = 3, nrow = length(Student_urn_size))
colnames(res_sys_urnings) = c("urn_size", "b_mse_sys", "HT")
res_sys_urnings[,1:2] = as.matrix(b_mse_sys_urnings)

res_urnings = matrix(0, ncol = 4, nrow = length(Student_urn_size) * length(Theta))
colnames(res_urnings) = c("urn_size", "ability", "b_mse", "HT")
res_urnings[,1:3] = as.matrix(b_mse_urnings[,1:3])

# run all K systems
for(nk in 1:length(Student_urn_size)){
  tic()
  #setting k
  us = Student_urn_size[nk]
  #getting baseline for k
  lq_baseline = b_mse_urnings[b_mse_urnings$urn_size == us, 4]
  uq_baseline =b_mse_urnings[b_mse_urnings$urn_size == us, 5]
  
  
  #status print
  print(paste0("Running urn size = ", us))
  
  #running the system
  res = urningsHT(n_students = n_students,
                  n_items = n_items,
                  n_reps = n_reps_baseline,
                  n_games = 1000,
                  student_urn_size = us,
                  item_urn_size = item_urn_size,
                  Theta = Theta,
                  Delta = Delta,
                  lq_baseline = lq_baseline,
                  uq_baseline = uq_baseline,
                  adaptive = 1,
                  paired = 1,
                  m_p = m_p,
                  s_p = sigma_p,
                  OS = "MAC")
  
  #saving the results to the output files 
  res_urnings[b_mse_urnings$urn_size == us, "HT"] = res$HT_student
  res_sys_urnings[nk, "HT"] = mean(res$HT_student)
  toc()
}

#saving the output
saveRDS(res_sys_urnings, "output/res2_sys_u3.RDS")
saveRDS(res_urnings, "output/res2_u3.RDS")



