#-------------------------------------------------------------------------------
#                         Baseline calculations for Urnings
# NOTE: if the abilties or the urn sizes that are used in the simulations change
# this matrix needs to be updated.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Dependencies
#-------------------------------------------------------------------------------
source("algo/baseline_functions.R")
source("sim2/sim2_settings.R")
source("util.R")

library(tidyverse)

Pi = expit(Theta_matrix[,4])
#-------------------------------------------------------------------------------
# containers
b_urnings = matrix(0, ncol = 7, nrow = length(Student_urn_size) * length(Theta))
colnames(b_urnings) = c("urn_size", "ability", "b_mse", "lq", "uq", "meanT", "varT")
b_sys_urnings = matrix(0, ncol = 2, nrow = length(Student_urn_size))
colnames(b_sys_urnings) = c("urn_size", "b_mse_sys")
#-------------------------------------------------------------------------------
# Calculation

for(us in 1:length(Student_urn_size)){
  curr_urn_size = Student_urn_size[us]
  print(curr_urn_size)
  res = baselineUrnings(Pi, curr_urn_size, n_reps = n_reps_baseline)
  b_sys_urnings[us,] = c(curr_urn_size,res$baseline_sys)
  b_urnings[((us-1)*length(Theta) + 1):(us*length(Theta)),] = res$baseline
}

#-------------------------------------------------------------------------------
# saving the output (RDS)
#-------------------------------------------------------------------------------
saveRDS(b_urnings, "output/b_u4.RDS")
saveRDS(b_sys_urnings, "output/b_sys_u4.RDS")



