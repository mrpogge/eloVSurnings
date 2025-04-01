################################################################################
# ELO algorithm wrapper:
# performs simulation with ELO algorithm 
################################################################################
eloHT = function(n_students,
                n_items,
                n_reps,
                n_games,
                K_s,
                K_i,
                Theta,
                Delta,
                lq_baseline,
                uq_baseline,
                adaptive=0,
                m_p=0,
                s_p=1,
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
  #   baseline: baseline MSE per student (numeric[n_students]) NOTE: this needs to be ordered as Theta
  #   baseline_sys: baseline MSE for the system (numeric)
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #
  # Output: 
  #   res: list[HT_sys, HT_student]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_HT.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_HT.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit!!
  
  theta_hat = rep(0,n_reps*n_students)
  delta_hat = rep(Delta,n_reps)
  mean_rating = rep(0,n_students*n_games)

  #-----------------------------------------------------------------------------
  # outcome objects
  HT_student = rep(n_games,n_students)
  is_HTv = rep(0,n_students)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("elo_HT",
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
           as.double(lq_baseline),
           as.double(uq_baseline),
           as.integer(adaptive),
           as.double(A),
           as.double(B),
           as.double(c(0:n_items+1)),
           as.integer(HT_student),
           as.integer(is_HTv),
           as.double(mean_rating))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(HT_sys = tmp[[17]] + 1)
  
}


################################################################################
# ELO algorithm wrapper:
# performs simulation with ELO algorithm 
################################################################################
eloDoubleHT = function(n_students,
                       n_items,
                       n_reps,
                       n_games,
                       K_s,
                       K_i,
                       Theta,
                       Delta,
                       lq_baseline,
                       uq_baseline,
                       m_p=0,
                       s_p=1,
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
  #   baseline: baseline MSE per student (numeric[n_students]) NOTE: this needs to be ordered as Theta
  #   baseline_sys: baseline MSE for the system (numeric)
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #
  # Output: 
  #   res: list[HT_sys, HT_student]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_HT.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_HT.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit!!
  
  theta_hat = rep(0,n_reps*n_students)
  delta_hat = rep(Delta,n_reps)
  mean_rating = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # outcome objects
  HT_student = rep(n_games,n_students)
  is_HTv = rep(0,n_students)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("elo_double_HT",
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
           as.double(lq_baseline),
           as.double(uq_baseline),
           as.double(A),
           as.double(B),
           as.double(c(0:n_items+1)),
           as.integer(HT_student),
           as.integer(is_HTv),
           as.double(mean_rating))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(HT_sys = tmp[[18]] + 1,
             mean_rating = matrix(tmp[[20]], nrow = n_students, ncol = n_games))
  
}


################################################################################
# ELO algorithm wrapper:
# performs simulation with ELO algorithm 
################################################################################
eloChange = function(n_students,
                     n_items,
                     n_reps,
                     n_games,
                     K_s,
                     K_i,
                     Theta,
                     Delta,
                     m_p=0,
                     s_p=1,
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
  #   baseline: baseline MSE per student (numeric[n_students]) NOTE: this needs to be ordered as Theta
  #   baseline_sys: baseline MSE for the system (numeric)
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #
  # Output: 
  #   res: list[HT_sys, HT_student]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_change.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_change.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit!!
  
  theta_hat = rep(0,n_reps*n_students)
  delta_hat = rep(Delta,n_reps)
  mean_rating = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # outcome objects
  bias_v = rep(0,n_students*n_games)
  var_v = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("elo_double_change",
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
           as.double(mean_rating),
           as.double(A),
           as.double(B),
           as.double(c(0:n_items+1)),
           as.double(bias_v),
           as.double(var_v))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(bias_v = matrix(tmp[[17]], nrow = n_students, ncol = n_games),
             var_v = matrix(tmp[[18]], nrow = n_students, ncol = n_games))
  
  if(OS == "MAC"){
    dyn.unload("algo/elo_change.so")
  } else if(OS == "WIN"){
    dyn.unload("algo/elo_change.dll")
  }
  
  return(res)
  
}

################################################################################
# ELO algorithm wrapper:
# performs simulation with ELO algorithm 
################################################################################
eloSimpleHT = function(n_students,
                       n_items,
                       n_reps,
                       n_games,
                       K_s,
                       Theta,
                       Delta,
                       lq_baseline,
                       uq_baseline,
                       m_p=0,
                       s_p=1,
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
  #   Theta: abilties (numeric[n_students))
  #   Delta: difficulties (numeric[n_items])
  #   baseline: baseline MSE per student (numeric[n_students]) NOTE: this needs to be ordered as Theta
  #   baseline_sys: baseline MSE for the system (numeric)
  #   adaptive: adaptivity (0-random or 1-adaptive)
  #   m_p: mean parameter for the Normal Kernel Method (numeric)
  #   s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  #   OS: operating system (MAC,WIN)
  #
  # Output: 
  #   res: list[HT_sys, HT_student]
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the compiled C routines
  #-------------------------------------------------------------------------------
  
  if(OS == "MAC"){
    dyn.load("algo/elo_HT.so")
  } else if(OS == "WIN"){
    dyn.load("algo/elo_HT.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  A = - 1/(2*s_p^2)
  B = m_p/s_p^2
  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit!!
  
  theta_hat = rep(0,n_reps*n_students)
  mean_rating = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # outcome objects
  HT_student = rep(n_games,n_students)
  is_HTv = rep(0,n_students)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("elo_simple_HT",
           as.double(K_s),
           as.integer(n_reps),
           as.integer(n_games),
           as.integer(n_students),
           as.integer(n_items),
           as.double(Theta),
           as.double(Delta),
           as.double(theta_hat),
           as.double(lq_baseline),
           as.double(uq_baseline),
           as.double(A),
           as.double(B),
           as.double(c(0:n_items+1)),
           as.integer(HT_student),
           as.integer(is_HTv),
           as.double(mean_rating))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(HT_sys = tmp[[14]] + 1,
             mean_rating = matrix(tmp[[16]], nrow = n_students, ncol = n_games))
  
}

