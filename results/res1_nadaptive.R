#-------------------------------------------------------------------------------
#                    Results of simulation 1
# Description: This file contains the results of the first simulation, including
# graphing and calculating summary statistics
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# dependencies
library(tidyverse)
library(plotly)
library(patchwork)
source("util.R")
#-------------------------------------------------------------------------------
# loading results

#elo
res_sys_elo = readRDS("output/res1_sys_elo.RDS")
res_elo = readRDS("output/res1_elo.RDS")

#urnings
res_sys_urnings = readRDS("output/res1_sys_urnings.RDS")
res_urnings = readRDS("output/res1_urnings.RDS")

#adding step size rank indicator to the results
res_sys_elo = res_sys_elo %>% as.data.frame() %>% mutate(step_size = 42-rank(K))
res_elo = res_elo %>% 
  as.data.frame() %>%
  left_join(res_sys_elo %>% select(K, step_size), by = "K") 
res_sys_urnings = res_sys_urnings %>% as.data.frame() %>% mutate(step_size = rank(urn_size))
res_urnings = res_urnings %>% 
  as.data.frame() %>%
  left_join(res_sys_urnings %>% select(urn_size, step_size), by = "urn_size")
#-------------------------------------------------------------------------------
# relationship of b_mse and step size (system level)

#create a single database
res_sys = matrix(0, ncol = 5, nrow = nrow(res_sys_elo) + nrow(res_sys_urnings))
colnames(res_sys) = c("system", "K_or_urn", "b_mse_sys", "HT", "step_size")
res_sys = as.data.frame(res_sys)

res_sys[1:nrow(res_sys_elo),1] = "elo"
res_sys[1:nrow(res_sys_elo),2:5] = res_sys_elo 
res_sys[(nrow(res_sys_elo) + 1):nrow(res_sys),1] = "urnings"
res_sys[(nrow(res_sys_elo) + 1):nrow(res_sys),2:5] = res_sys_urnings

# plotting the results

res_sys %>%
  ggplot(aes(x = b_mse_sys, y = HT, color = system)) +
  geom_line() +
  labs(x = "MSE",
       y = "HT") +
  theme_minimal()

#-------------------------------------------------------------------------------
# matched systems 

#calculate distance between baseline mse-s
distances = matrix(0, ncol = 4, nrow = nrow(res_sys_elo) * nrow(res_sys_urnings))
colnames(distances) = c("K", "urn_size", "b_mse_diff", "HT_diff")
counter = 1
for(i in 1:nrow(res_sys_elo)){
  for(j in 1:nrow(res_sys_urnings)){
    distances[counter,1] = res_sys_elo$K[i]
    distances[counter,2] = res_sys_urnings$urn_size[j]
    distances[counter,3] = abs(res_sys_elo$b_mse_sys[i] - res_sys_urnings$b_mse_sys[j])
    distances[counter,4] = res_sys_urnings$HT[j] - res_sys_elo$HT[i]
    counter = counter + 1
  }
}
distances = distances %>% as.data.frame()

matched_pairs = matrix(0, ncol = 4, nrow = nrow(res_sys_elo))
colnames(matched_pairs) = c("urn_size", "matched_K", "matched_HT", "b_mse_diff")
for(i in 1:nrow(matched_pairs)){
  us = res_sys_urnings$urn_size[i]
  min_bmse = min(distances[distances[,"urn_size"] == us, "b_mse_diff"])
  min_K_HT = distances[distances[,"urn_size"] == us & distances[,"b_mse_diff"] == min_bmse, c("K", "HT_diff")]
  matched_pairs[i,1] = as.integer(us)
  matched_pairs[i,2:3] = as.numeric(min_K_HT)
  matched_pairs[i,4] = as.numeric(min_bmse)
}
matched_pairs = matched_pairs %>% as.data.frame()
  
#plot urn size vs matched HT with a secondary axes for matched K
p1 = matched_pairs %>%
  ggplot(aes(x = urn_size, y = matched_HT, color = matched_K)) + 
  geom_point() + 
  labs(x = "urn size",
       y = "HT difference") +
  theme_minimal()

#-------------------------------------------------------------------------------
# ht contour plot for studentwise ht-s
#elo
res_elo_logit1 = res_elo %>%
  filter(ability < 2 & ability > 1) %>%
  group_by(K,step_size) %>%
  summarise(HT = mean(HT),
            mMSE = mean(b_mse))

res_urnings_logit1 = res_urnings %>%
  filter(logit(ability) < 2 & logit(ability) > 1) %>%
  group_by(urn_size, step_size) %>%
  summarise(HT = mean(HT),
            mMSE = mean(b_mse))

colnames(res_urnings_logit1)[1] = "K"
res_logit1 = cbind(rep(c("elo", "urnings"), each = nrow(res_sys_elo)), rbind(res_elo_logit1, res_urnings_logit1))
colnames(res_logit1)[1] = "system"
res_logit1 %>%
  as.data.frame() %>%
  ggplot(aes(x = mMSE, y = HT, color = system)) +
  geom_point() +
  geom_line()
