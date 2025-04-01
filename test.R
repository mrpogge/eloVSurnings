#calculate theoretical MSE of Urnings with a given ability distribution and 
# urn size

#################################################################################
# fix values 
n_students = 900
n_items = 200
n_games = 500
student_urn_sizes = seq(2,512,by=4)

# generate theta as equal quantiles of a normal distribution with 0 mean and 1 sd
Theta = qnorm(seq(1/(n_students+1),n_students/(n_students+1),length=n_students),0,1)
sample(Theta, n_students)
Delta = qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1.5)

# calculate the variances of the Urnings theoretical distribution
Pi = exp(Theta) / (1 + exp(Theta))

########################################################################################################
# running 1 elo system test
#######################################################################################################

source("algo/elo_wrapper.R")
#-----------------#
# setting basic parameters 

n = 900
m = 200
games = 500
reps = 10

K_s = seq(0.03,1, by = 0.02)
mse_elo = matrix(0, nrow = length(K_s), ncol = 2)
mse_elo[,1] = K_s
colnames(mse_elo) = c("K", "meanMSE")
#-----------------#
# running elo function
counter = 1
for(k in K_s){
  print(k)
  a = elo(n = n,
          m = m,
          games = games,
          reps = reps,
          fixed_items = 0,
          adaptive = 0,
          K = k,
          theta = Theta)
  
  mse_elo[counter, 2] = mean(colMeans(a$mse_theta)[(games-99):games])
  counter = counter + 1
}
library(ggplot2)
mse_elo %>%
  as.data.frame() %>%
  ggplot(aes(x = K, y = meanMSE)) +
  geom_line() +
  labs(title = "Elo system: K vs MSE",
       x = "K",
       y = "Mean MSE") +
  theme_minimal()


#################################################################################
# urn size and mse relationship 
urn_sizes = seq(5, 100, by = 1)
mse_s = matrix(0, nrow = length(Pi) * length(urn_sizes), ncol = 3)

counter = 1
for(i in 1:length(Pi)){
  for(j in 1:length(urn_sizes)){
    mse_s[counter, ] = c(Pi[i], urn_sizes[j], 1/urn_sizes[j] *Pi[i] * (1 - Pi[i]))
    counter = counter + 1
  }
}
colnames(mse_s) = c("theta", "urn_size", "mse")

library(tidyverse)

mse_s %>%
  as.data.frame() %>%
  group_by(urn_size) %>%
  summarise(meanMSE = mean(mse)) %>%
  ggplot(aes(x = urn_size, y = meanMSE)) +
  geom_line() 

mse_urnings = mse_s %>%
  as.data.frame() %>%
  group_by(urn_size) %>%
  filter(urn_size %in% c(8,16,32,64)) %>%
  summarise(meanMSE = mean(mse))

#find k where the mse_urnings and the mse_elo difference is minimal given an urn size
differences = outer(mse_urnings$meanMSE,mse_elo[,"meanMSE"], "-") %>%
  abs() 

colnames(differences) = mse_elo[,"K"]
rownames(differences) = mse_urnings$urn_size

diff_long = matrix(0,ncol = 3, nrow = nrow(differences) * ncol(differences))
colnames(diff_long) = c("urn_size", "K", "mse")
for(i in 1:length(mse_urnings$urn_size)){
  for(j in 1:length(mse_elo[,"K"])){
    diff_long[(i-1) * length(mse_elo[,"K"]) + j,] = c(mse_urnings$urn_size[i], mse_elo[,"K"][j], differences[i,j])
  }
}

diff_long %>%
  as.data.frame() %>%
  ggplot(aes(x = urn_size, y = K, fill = mse)) +
  geom_tile() +
  labs(title = "Urn size vs K",
       x = "Urn size",
       y = "K",
       fill = "Difference") +
  theme_minimal()

selections = matrix(0, nrow = nrow(differences), ncol = 3)
for(r in 1:nrow(differences)){
  print(paste("urn size: ", rownames(differences)[r], " K: ", colnames(differences)[which.min(differences[r,])]))
  step_size = as.numeric(c(rownames(differences)[r], colnames(differences)[which.min(differences[r,])]))
  selections[r,1:2] = step_size
  selections[r,3] = differences[r,which.min(differences[r,])]
}
colnames(selections) = c("urn_size", "K", "mse")

selections %>%
  as.data.frame() %>%
  ggplot(aes(x = as.factor(urn_size), y = K)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0.1,1,by=0.1)) +
  labs(title = "",
       x = "Urn size",
       y = "K") +
  jtools::theme_apa(legend.font.size = 10)

selections %>%
  as.data.frame() %>%
  ggplot(aes(x = urn_size, y = mse)) +
  geom_point() +
  labs(title = "",
       x = "Urn size",
       y = "MSE") +
  jtools::theme_apa(legend.font.size = 10)

#check the precalced baselines
b_mse_sys_elo = readRDS("output/b_mse_sys_elo.RDS")
b_mse_sys_urnings = readRDS("output/b_mse_sys_urnings.RDS")

b_mse_sys_elo %>%
  as.data.frame() %>%
  ggplot(aes(x = max(K)-K, y = b_mse_sys)) +
  geom_line() +
  labs(title = "Elo system: K vs MSE",
       x = "K",
       y = "Mean MSE") +
  theme_minimal() 

b_mse_sys_urnings %>%
  as.data.frame() %>%
  ggplot(aes(x = urn_size, y = b_mse_sys)) +
  geom_line() +
  labs(title = "Urn size: urn size vs MSE",
       x = "Urn size",
       y = "Mean MSE") +
  theme_minimal()

#compare the two systems
differences = outer(b_mse_sys_urnings$b_mse_sys,b_mse_sys_elo$b_mse_sys, "-") %>%
  abs()

colnames(differences) = b_mse_sys_elo$K
rownames(differences) = b_mse_sys_urnings$urn_size

diff_long = matrix(0,ncol = 3, nrow = nrow(differences) * ncol(differences))
colnames(diff_long) = c("urn_size", "K", "mse")
for(i in 1:length(b_mse_sys_urnings$urn_size)){
  for(j in 1:length(b_mse_sys_elo$K)){
    diff_long[(i-1) * length(b_mse_sys_elo$K) + j,] = c(b_mse_sys_urnings$urn_size[i], b_mse_sys_elo$K[j], differences[i,j])
  }
}

selections = matrix(0, nrow = nrow(differences), ncol = 3)
for(r in 1:nrow(differences)){
  print(paste("urn size: ", rownames(differences)[r], " K: ", colnames(differences)[which.min(differences[r,])]))
  step_size = as.numeric(c(rownames(differences)[r], colnames(differences)[which.min(differences[r,])]))
  selections[r,1:2] = step_size
  selections[r,3] = differences[r,which.min(differences[r,])]
}
colnames(selections) = c("urn_size", "K", "mse")

selections %>%
  as.data.frame() %>%
  ggplot(aes(x = urn_size, y = K)) +
  geom_line() +
  labs(title = "Urn size vs K",
       x = "Urn size",
       y = "K") +
  theme_minimal()

#--------------------------------------------------------------------------------
# test elo HT function

#get a baseline
b_mse_sys_elo = readRDS("output/b_sys_elo.RDS")
b_mse_elo = readRDS("output/b_elo.RDS")

K = b_mse_sys_elo[3,1]
K_i = K_item
lq_baseline = b_mse_elo[b_mse_elo$K_s == K, 4]
uq_baseline = b_mse_elo[b_mse_elo$K_s == K, 5]
K_s = K
adaptive=0
m_p=0
s_p=1
OS = "MAC"
n_reps = 1000
n_games = 1000
#--------------------------------------------------------------------------------
# test urnings HT function

#get a baseline
b_mse_sys_urnings = readRDS("output/b_sys_urnings.RDS")
b_mse_urnings = as.data.frame(readRDS("output/b_urnings.RDS"))

student_urn_size = 32
lq_baseline = b_mse_urnings[b_mse_urnings$urn_size == student_urn_size, 4]
uq_baseline =b_mse_urnings[b_mse_urnings$urn_size == student_urn_size, 5]
adaptive = 0
m_p = 0
s_p = 1
OS = "MAC"
n_reps = 1000
n_games = 1000

#-------------------------------------------------------------------------------
# CLT approximation to find a boundary

sampling_for_mean = numeric(length = 1000)
for(r in 1:n_reps){
sampling_for_mean[r] = mean(rbinom(1000, 32, 0.5))
}

mean(sampling_for_mean)
var(sampling_for_mean)

r1 = rbinom(1000, 32, 0.5)
mean(r1)
var(r1) 

true_mean = 0.5*32
true_var = 0.5*0.5*32

#get the 95% CI for the mean by the normal approximation
df = length(sampling_for_mean) - 1
alpha = 0.05
t.score = qt(p=alpha/2, df=df,lower.tail=F)
se = sqrt(true_var/length(sampling_for_mean))
CI = c(true_mean - t.score * se, true_mean + t.score * se)

is_converged = numeric(length(nrow(a)))
for(i in 1:nrow(a)){
  check1 = mean(a[i,])
  theta_check1 = expit(Theta[i]) * 32
  var_check1 = (expit(Theta[i]) * (1-expit(Theta[i])))*32
  
  df = 999
  alpha = 0.05
  t.score = qt(p=alpha/2, df=df,lower.tail=F)
  se = sqrt(var_check1/1000)
  CI = c(theta_check1 - t.score * se, theta_check1 + t.score * se)
  is_converged[i] = as.numeric(check1>CI[1] & check1<CI[2])
}


#--------------------------------------------------------------------------------
# change test
student_urn_size = 64
adaptive = 1
m_p = 0
s_p = 1
OS = "MAC"
n_reps = 100
n_games = 500
Theta = matrix(0, nrow = n_students, ncol = n_games)
Theta[,1] =  qnorm(seq(1/(n_students+1),n_students/(n_students+1),length=n_students),mu_theta,sigma_theta) 
for(i in 2:n_games){
  Theta[,i] = Theta[,i-1] + 1/500
}
paired = 1


#--------------------------------------------------------------------------------
# change test
K_s = 0.4
adaptive = 1
m_p = 0
s_p = 1
OS = "MAC"
n_reps = 100
n_games = 10
Theta = matrix(0, nrow = n_students, ncol = n_games)
Theta[,1] =  qnorm(seq(1/(n_students+1),n_students/(n_students+1),length=n_students),mu_theta,sigma_theta) 
for(i in 2:n_games){
  Theta[,i] = Theta[,i-1] + 1/n_games
}
paired = 1
K_i = 0.3

res = rbind(res_elo, res_urnings)
res = as.data.frame(res)
res = cbind(res,rep(c("ers", "urnings"), each = 41))
colnames(res)[5] = "algo"

library(tidyverse)

res %>%
  ggplot(aes(x = var, y = bias, color = algo)) +
  geom_point() +
  labs(x = "Urn size",
       y = "Mean MSE") +
  theme_minimal()