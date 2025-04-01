#------------------------------------------------------------------------------
# Settings for the 3rd simulation 
#------------------------------------------------------------------------------
set.seed(13181917)
# Dependencies
source("util.R")
#-------------------------------------------------------------------------------
# General system settings

n_students = 1000
n_items = 200
n_games = 1000
n_reps = 200 #1000 in the real simulation

#-------------------------------------------------------------------------------
# General baseline settings
n_games_baseline = 400
n_reps_baseline = 200 #2000 in the real simulation


Mu_theta = c(0, 0, logit(0.7), logit(0.9))
sigma_theta = 1.5 

Mu_delta = 0
sigma_delta = 2

Theta_matrix = matrix(0, n_students, length(Mu_theta))
for(i in 1:length(Mu_theta)){
  theta_middle = qnorm(seq(1/(n_students+1),n_students/(n_students+1),length=n_students),Mu_theta[i],sigma_theta)
  Theta_matrix[,i] = sample(theta_middle, length(theta_middle), replace = FALSE)
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
window_size = 200
