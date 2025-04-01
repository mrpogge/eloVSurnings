#-------------------------------------------------------------------------------
#                           Urnings wrapper functions
# Description: This file contains functions which wrap the simulation code for Urnings
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#dependencies

source("util.R")

urningsHT = function(n_students,
                     n_items,
                     n_reps,
                     n_games,
                     student_urn_size,
                     item_urn_size,
                     Theta,
                     Delta,
                     lq_baseline,
                     uq_baseline,
                     adaptive=0,
                     m_p=0,
                     s_p=1,
                     paired = 0,
                     OS = "MAC"){
#--------------------------------------------------------------------------------
# arguments
# n_students: number of students (int)
# n_items: number of items (int)
# n_reps: number of replications (int)
# n_games: number of games (int)
# student_urn_size: urn size for students (int)
# item_urn_size: urn size for items (int)
# Theta: abilties (numeric[n_students))
# Delta: difficulties (numeric[n_items])
# baseline: baseline MSE per student (numeric[n_students]) NOTE: this needs to be ordered as Theta
# baseline_sys: baseline MSE for the system (numeric)
# adaptive: adaptivity (0-random or 1-adaptive)
# m_p: mean parameter for the Normal Kernel Method (numeric)
# s_p: standard deviation parameter for the Normal Kernel Method (numeric)
# OS: operating system (MAC,WIN)
  
#return
# res: list[HT_sys, HT_student]
#-------------------------------------------------------------------------------

  #-------------------------------------------------------------------------------
  # loading the C files
  
  if(OS == "MAC"){
    dyn.load("algo/urnings_HT.so")
  } else if(OS == "WIN"){
    dyn.load("algo/urnings_HT.dll")
  }

  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  
  Prob2=matrix(1,nrow=student_urn_size+1,ncol=item_urn_size+1)
  if(adaptive == 1){ 
    for(i in 0:(student_urn_size)){
      for(j in 0:(item_urn_size)){
        l=log((i+1)/(student_urn_size-i+1))-log((j+1)/(item_urn_size-j+1))
        Prob2[i+1,j+1]=dnorm(1/(1+exp(-l)),m_p,s_p)
      }
    }  
  }
  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit (0.5 probability)
  
  theta_hat = rep(as.integer(0.5*student_urn_size),n_reps*n_students)
  delta_hat = rep(rbinom(n_items, item_urn_size, expit(Delta)), n_reps)
  mean_rating = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # aux for paired update
  queue = rep(0, n_items*n_reps)
  P_queue = rep(0,n_reps)
  
  #-----------------------------------------------------------------------------
  # outcome objects
  HT_student = rep(n_games,n_students)
  is_HTv = rep(0,n_students)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("urnings_HT",
           as.integer(student_urn_size),
           as.integer(item_urn_size),
           as.integer(n_reps),
           as.integer(n_games),
           as.integer(n_students),
           as.integer(n_items),
           as.double(Theta),
           as.double(Delta),
           as.integer(theta_hat),
           as.integer(delta_hat),
           as.double(lq_baseline),
           as.double(uq_baseline),
           as.integer(adaptive),
           as.double(Prob2),
           as.double(c(0:n_items+1)),
           as.integer(paired),
           as.integer(queue),
           as.integer(P_queue),
           as.integer(HT_student),
           as.integer(is_HTv),
           as.double(mean_rating))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(HT_student = tmp[[19]] + 1,
             M_i = matrix(tmp[[21]], nrow = n_students, ncol = n_games))
  
  return(res)
}

urningsSimpleHT = function(n_students,
                     n_items,
                     n_reps,
                     n_games,
                     student_urn_size,
                     Theta,
                     Delta,
                     lq_baseline,
                     uq_baseline,
                     m_p=0,
                     s_p=1,
                     OS = "MAC"){
  #--------------------------------------------------------------------------------
  # arguments
  # n_students: number of students (int)
  # n_items: number of items (int)
  # n_reps: number of replications (int)
  # n_games: number of games (int)
  # student_urn_size: urn size for students (int)
  # item_urn_size: urn size for items (int)
  # Theta: abilties (numeric[n_students))
  # Delta: difficulties (numeric[n_items])
  # baseline: baseline MSE per student (numeric[n_students]) NOTE: this needs to be ordered as Theta
  # baseline_sys: baseline MSE for the system (numeric)
  # adaptive: adaptivity (0-random or 1-adaptive)
  # m_p: mean parameter for the Normal Kernel Method (numeric)
  # s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  # OS: operating system (MAC,WIN)
  
  #return
  # res: list[HT_sys, HT_student]
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the C files
  
  if(OS == "MAC"){
    dyn.load("algo/urnings_HT.so")
  } else if(OS == "WIN"){
    dyn.load("algo/urnings_HT.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  
  Prob2=matrix(1,nrow=student_urn_size+1,ncol=n_items)
  for(i in 0:(student_urn_size)){
    for(j in 1:n_items){
      l=log((i+1)/(student_urn_size-i+1))-Delta[j]
      Prob2[i+1,j]=dnorm(1/(1+exp(-l)),m_p,s_p)
    }
  }  

  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit (0.5 probability)
  
  theta_hat = rep(as.integer(0.5*student_urn_size),n_reps*n_students)
  mean_rating = rep(0,n_students*n_games)
  
  
  #-----------------------------------------------------------------------------
  # outcome objects
  HT_student = rep(n_games,n_students)
  is_HTv = rep(0,n_students)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("urnings_simple_HT",
           as.integer(student_urn_size),
           as.integer(n_reps),
           as.integer(n_games),
           as.integer(n_students),
           as.integer(n_items),
           as.double(Theta),
           as.double(Delta),
           as.integer(theta_hat),
           as.double(lq_baseline),
           as.double(uq_baseline),
           as.double(Prob2),
           as.double(c(0:n_items+1)),
           as.integer(HT_student),
           as.integer(is_HTv),
           as.double(mean_rating))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(HT_student = tmp[[13]] + 1,
             M_i = matrix(tmp[[15]], nrow = n_students, ncol = n_games))
  
  return(res)
}


urningsChange = function(n_students,
                         n_items,
                         n_reps,
                         n_games,
                         student_urn_size,
                         item_urn_size,
                         Theta,
                         Delta,
                         adaptive=0,
                         m_p=0,
                         s_p=1,
                         paired = 0,
                         OS = "MAC"){
  #--------------------------------------------------------------------------------
  # arguments
  # n_students: number of students (int)
  # n_items: number of items (int)
  # n_reps: number of replications (int)
  # n_games: number of games (int)
  # student_urn_size: urn size for students (int)
  # item_urn_size: urn size for items (int)
  # Theta: abilties (numeric[n_students))
  # Delta: difficulties (numeric[n_items])
  # adaptive: adaptivity (0-random or 1-adaptive)
  # m_p: mean parameter for the Normal Kernel Method (numeric)
  # s_p: standard deviation parameter for the Normal Kernel Method (numeric)
  # OS: operating system (MAC,WIN)
  
  #return
  # res: list[HT_sys, HT_student]
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # loading the C files
  
  if(OS == "MAC"){
    dyn.load("algo/urnings_change.so")
  } else if(OS == "WIN"){
    dyn.load("algo/urnings_change.dll")
  }
  
  #-----------------------------------------------------------------------------
  # setting adaptive item selection params
  
  Prob2=matrix(1,nrow=student_urn_size+1,ncol=item_urn_size+1)
  if(adaptive == 1){ 
    for(i in 0:(student_urn_size)){
      for(j in 0:(item_urn_size)){
        l=log((i+1)/(student_urn_size-i+1))-log((j+1)/(item_urn_size-j+1))
        Prob2[i+1,j+1]=dnorm(1/(1+exp(-l)),m_p,s_p)
      }
    }  
  }
  
  #-----------------------------------------------------------------------------
  # starting values  and containers: NOTE: start at 0 logit (0.5 probability)
  
  theta_hat = rep(as.integer(0.5*student_urn_size),n_reps*n_students)
  delta_hat = rep(rbinom(n_items, item_urn_size, expit(Delta)), n_reps)
  mean_rating = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # aux for paired update
  queue = rep(0, n_items*n_reps)
  P_queue = rep(0,n_reps)
  
  #-----------------------------------------------------------------------------
  # outcome objects
  bias_v = rep(0,n_students*n_games)
  var_v = rep(0,n_students*n_games)
  
  #-----------------------------------------------------------------------------
  # run the simulation
  tmp = .C("urnings_change",
           as.integer(student_urn_size),
           as.integer(item_urn_size),
           as.integer(n_reps),
           as.integer(n_games),
           as.integer(n_students),
           as.integer(n_items),
           as.double(Theta),
           as.double(Delta),
           as.integer(theta_hat),
           as.integer(delta_hat),
           as.integer(adaptive),
           as.double(Prob2),
           as.double(c(0:n_items+1)),
           as.integer(paired),
           as.integer(queue),
           as.integer(P_queue),
           as.double(mean_rating),
           as.double(bias_v),
           as.double(var_v))
  
  #-----------------------------------------------------------------------------
  # output
  res = list(bias_v = matrix(tmp[[18]], nrow = n_students, ncol = n_games),
             var_v = matrix(tmp[[19]], nrow = n_students, ncol = n_games))
  
  if(OS == "MAC"){
    dyn.unload("algo/urnings_change.so")
  } else if(OS == "WIN"){
    dyn.unload("algo/urnings_change.dll")
  }
  
  return(res)
}