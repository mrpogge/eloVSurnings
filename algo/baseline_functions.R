#-------------------------------------------------------------------------------
#                           Baseline MSE
# Description: This file contains functions to calculate the baseline MSE of Urnings
# and Elo with different K/urn_size and abilities.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Baseline simulation function for Urnings
#-------------------------------------------------------------------------------

baselineUrnings = function(Pi,urn_size, n_reps = 1000){
  #-----------------------------------------------------------------------------
  # baselineUrnings_p: calculate the baseline MSE (95SI) of Urnings based on the urn size 
  # and the expit transformed ability for a single person. And the 95%SI of the limiting distribtion
  # for each ability
  #
  # Arguments: 
  #   pi: expit transformed ability of a person (numeric or numeric[length(pi)])
  #   urn_size: number of items in the urn (int or int[length(pi)]
  #
  # Output: 
  #   MSE of Urnings (numeric or numeric[length(pi)]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # sampling from the limiting distribution
  n_students = length(Pi)
  container = matrix(0, nrow = n_students, ncol = n_reps)
  for(r in 1:n_reps){
    u = rbinom(n_students, urn_size, Pi)
    container[,r] = u
  }
  
  #-----------------------------------------------------------------------------
  # mse of each student's limiting distribution
  mse = rowMeans(((container / urn_size) - Pi)^2)
  theoretical_mse= mean(Pi * (1-Pi) * 1/urn_size)
  #-----------------------------------------------------------------------------
  # 95CI for mean ratings
  th_mean = Pi * urn_size
  th_var = Pi * (1-Pi) * urn_size
  lq = th_mean + qnorm(0.025) * sqrt(th_var/n_reps)
  uq = th_mean + qnorm(0.975) * sqrt(th_var/n_reps)
  #-----------------------------------------------------------------------------
  # outcome
  baseline = matrix(0, nrow = n_students, ncol = 7)
  colnames(baseline) = c("urn_size", "ability", "mse", "lq", "uq", "meanT", "varT")
  baseline[,1] = rep(urn_size, n_students)
  baseline[,2] = Pi
  baseline[,3] = mse
  baseline[,4] = lq
  baseline[,5] = uq
  baseline[,6] = th_mean
  baseline[,7] = th_var
  
  return(list(baseline_sys = theoretical_mse, baseline = baseline))
}

#-------------------------------------------------------------------------------
# Baseline simulation function for Elo
#-------------------------------------------------------------------------------
baselineElo = function(n_students,
                       n_items,
                       n_reps,
                       n_games,
                       K_s,
                       K_i,
                       Theta,
                       Delta,
                       adaptive = 0,
                       m_p = 0,
                       s_p = 1,
                       OS = "MAC"){
  #-----------------------------------------------------------------------------
  # baselineElo: calculate the mean and variance of MSE of Elo based on the K and
  # the abilties. NOTE: this require some further analysis to decide where MSE is stable.
  #
  # Arguments:
  #   n_students: number of students (int)
  #   n_items: number of items (int)
  #   n_reps: number of replications (int)
  #   n_games: number of games (int)
  #   K_s: K for students (numeric)
  #   K_i: K for items (numeric)
  #   Theta: abilties (numeric[n_students))
  #   Delta: difficulties (numeric[n_items])
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #   seeds: seeds for replications (numeric[n_reps])
  #
  # Output: 
  #   res: [meanMSE matrix[n_students, n_games]]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_baseline.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_baseline.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers
  #NOTE: we start the algorithm for their true values so the algorithm converges faster
  
  theta_hat = rep(Theta,n_reps)
  delta_hat = rep(Delta,n_reps)
  
  mean_mse = rep(0,n_students*n_games)
  
  tmp<-.C("elo_baseline",
          as.double(K_s),
          as.double(K_i),
          as.integer(n_reps),
          as.integer(n_games),
          as.integer(n_students),
          as.integer(n_items),
          as.double(Theta),
          as.double(Delta),
          as.double(theta_hat),
          as.double(delta_hat),
          as.double(mean_mse),
          as.double(mean_mse),
          as.double(mean_mse),
          as.integer(adaptive),
          as.double(A),
          as.double(B),
          as.double(c(0:n_items+1)))
  
  res = list(Theta = Theta, 
             MSE_i = matrix(tmp[[11]],nrow = n_students, ncol = n_games),
             meanT = matrix(tmp[[12]],nrow = n_students, ncol = n_games),
             varT = matrix(tmp[[13]],nrow = n_students, ncol = n_games))
  return(res)
}

#-------------------------------------------------------------------------------
# Baseline simulation function for double Elo
#-------------------------------------------------------------------------------
baselineDoubleElo = function(n_students,
                             n_items,
                             n_reps,
                             n_games,
                             K_s,
                             K_i,
                             Theta,
                             Delta,
                             m_p = 0,
                             s_p = 1,
                             OS = "MAC"){
  #-----------------------------------------------------------------------------
  # baselineDoubleElo: calculate the mean and variance of MSE of Elo based on the K and
  # the abilties. NOTE: this require some further analysis to decide where MSE is stable.
  #
  # Arguments:
  #   n_students: number of students (int)
  #   n_items: number of items (int)
  #   n_reps: number of replications (int)
  #   n_games: number of games (int)
  #   K_s: K for students (numeric)
  #   K_i: K for items (numeric)
  #   Theta: abilties (numeric[n_students))
  #   Delta: difficulties (numeric[n_items])
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #   seeds: seeds for replications (numeric[n_reps])
  #
  # Output: 
  #   res: [meanMSE matrix[n_students, n_games]]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_baseline.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_baseline.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers
  #NOTE: we start the algorithm for their true values so the algorithm converges faster
  
  theta_hat = rep(Theta,n_reps)
  delta_hat = rep(Delta,n_reps)
  
  mean_mse = rep(0,n_students*n_games)
  
  tmp<-.C("elo_double_baseline",
          as.double(K_s),
          as.double(K_i),
          as.integer(n_reps),
          as.integer(n_games),
          as.integer(n_students),
          as.integer(n_items),
          as.double(Theta),
          as.double(Delta),
          as.double(theta_hat),
          as.double(delta_hat),
          as.double(theta_hat),
          as.double(delta_hat),
          as.double(mean_mse),
          as.double(mean_mse),
          as.double(mean_mse),
          as.double(A),
          as.double(B),
          as.double(c(0:n_items+1)))
  
  res = list(Theta = Theta, 
             MSE_i = matrix(tmp[[13]],nrow = n_students, ncol = n_games),
             M_i = matrix(tmp[[14]],nrow = n_students, ncol = n_games),
             Var_i = matrix(tmp[[15]],nrow = n_students, ncol = n_games))
  return(res)
}

#-------------------------------------------------------------------------------
# Baseline simulation function for simple Elo
#-------------------------------------------------------------------------------
baselineSimpleElo = function(n_students,
                             n_items,
                             n_reps,
                             n_games,
                             K_s,
                             Theta,
                             Delta,
                             m_p = 0,
                             s_p = 1,
                             OS = "MAC"){
  #-----------------------------------------------------------------------------
  # baselineDoubleElo: calculate the mean and variance of MSE of Elo based on the K and
  # the abilties. NOTE: this require some further analysis to decide where MSE is stable.
  #
  # Arguments:
  #   n_students: number of students (int)
  #   n_items: number of items (int)
  #   n_reps: number of replications (int)
  #   n_games: number of games (int)
  #   K_s: K for students (numeric)
  #   K_i: K for items (numeric)
  #   Theta: abilties (numeric[n_students))
  #   Delta: difficulties (numeric[n_items])
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #   seeds: seeds for replications (numeric[n_reps])
  #
  # Output: 
  #   res: [meanMSE matrix[n_students, n_games]]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_baseline.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_baseline.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers
  #NOTE: we start the algorithm for their true values so the algorithm converges faster
  
  theta_hat = rep(Theta,n_reps)
  
  mean_mse = rep(0,n_students*n_games)
  
  tmp<-.C("elo_simple_baseline",
          as.double(K_s),
          as.integer(n_reps),
          as.integer(n_games),
          as.integer(n_students),
          as.integer(n_items),
          as.double(Theta),
          as.double(Delta),
          as.double(theta_hat),
          as.double(mean_mse),
          as.double(mean_mse),
          as.double(mean_mse),
          as.double(A),
          as.double(B),
          as.double(c(0:n_items+1)))
  
  res = list(Theta = Theta, 
             MSE_i = matrix(tmp[[9]],nrow = n_students, ncol = n_games),
             M_i = matrix(tmp[[10]],nrow = n_students, ncol = n_games),
             Var_i = matrix(tmp[[11]],nrow = n_students, ncol = n_games))
  return(res)
}
