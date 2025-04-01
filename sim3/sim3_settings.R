#------------------------------------------------------------------------------
# Settings for the 3rd simulation 
#------------------------------------------------------------------------------

# Dependencies
source("util.R")
#-------------------------------------------------------------------------------
# General system settings

n_students = 1000
n_items = 200
n_games = 501
n_reps = 200 #1000 in the real simulation

#-------------------------------------------------------------------------------
# population settings

Mu_theta = c(0, 0, logit(0.7), logit(0.9))
sigma_theta = 1.5 

Mu_delta = 0
sigma_delta = 2

#change settings
change_amounts = rep(c(2,2,2,2,2),times = n_students / 5) / (n_games-1) #student change amount allocation equally

#latent change is set s.t. the 251th iteration's distribution matches the normal with the given settings
Theta_array = array(0, dim = c(n_students, n_games,length(Mu_theta)))
for(i in 1:length(Mu_theta)){
  
  theta_middle = qnorm(seq(1/(n_students+1),n_students/(n_students+1),length=n_students),Mu_theta[i],sigma_theta)
  Theta_array[,251,i] = sample(theta_middle, length(theta_middle), replace = FALSE)
  
  for(j in 1:((n_games-1)/2)){
    k = j + 251
    l = 251 -j
    Theta_array[,l,i] = Theta_array[,l+1,i] - change_amounts
    Theta_array[,k,i] = Theta_array[,k-1,i] + change_amounts
  }
}

#item difficulties
Delta = qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),Mu_delta,sigma_delta)

#-------------------------------------------------------------------------------
# starting values are noninformed, but could be changed later 

#TODO!!

#-------------------------------------------------------------------------------
# adaptivity settings

Mu_p = c(0, logit(0.7), logit(0.7), logit(0.9))
sigma_p = 1

#-------------------------------------------------------------------------------
#student urn sizes
Student_urn_size = seq(10, 90, by = 2)

#student K values
K_s = seq(0.1,0.8, length.out = 41)

#TODO: check the validity of such selection, for equal item variances.
item_urn_size = 45
K_i = 0.2

#-------------------------------------------------------------------------------
# util

# window size for the baseline calculation
window_size = 200

