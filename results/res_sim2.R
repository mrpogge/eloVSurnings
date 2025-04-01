#-------------------------------------------------------------------------------
# Results of the second simulation
#-------------------------------------------------------------------------------

#dependencies
library(tidyverse)
source("util.R")

#-------------------------------------------------------------------------------
# loading the system level results of the simulation

#elo
b_elo1 = readRDS("output/res2_sys_e1.RDS")
b_elo2 = readRDS("output/res2_sys_e2.RDS")
b_elo3 = readRDS("output/res2_sys_e3.RDS")
b_elo4 = readRDS("output/res2_sys_e4.RDS")

#urnings
b_urnings1 = readRDS("output/res2_sys_u1.RDS")
b_urnings2 = readRDS("output/res2_sys_u2.RDS")
b_urnings3 = readRDS("output/res2_sys_u3.RDS")
b_urnings4 = readRDS("output/res2_sys_u4.RDS")

# adding case indicator variables
algo = rep(c("elo", "urnings"), each = 164)
case = rep(c(1,2,3,4,1,2,3,4), each = 41)

dat = rbind(b_elo1, b_elo2, b_elo3, b_elo4, b_urnings1, b_urnings2, b_urnings3, b_urnings4) 
dat = as.data.frame(dat)
dat = cbind(dat, algo, case)

#-------------------------------------------------------------------------------
#visualise results:
# MSE vs HT for each case facet wrapped add a legend

plt1 = dat %>%
  ggplot(aes(x = b_mse_sys, y = HT, color = algo)) +
  geom_point() +
  facet_wrap(~case) +
  labs(x = "MSE",
       y = "HT") +
  theme_minimal()
