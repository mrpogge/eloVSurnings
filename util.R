#-------------------------------------------------------------------------------
#                         Utility functions
#-------------------------------------------------------------------------------

expit = function(x){
  #-----------------------------------------------------------------------------
  # expit: expit transformation (logit -> "probability")
  # 
  # Arguments:
  #   x: numeric (or vector)
  #
  # Output:
  #   numeric (or vector)
  return(1/(1 + exp(-x)))
}

logit = function(x){
  #-----------------------------------------------------------------------------
  # logit: logit transformation ("probability" -> logit)
  # 
  # Arguments:
  #   x: numeric (or vector)
  #
  # Output:
  #   numeric (or vector)
  return(log(x/(1-x)))
}

seq_stationarity = function(x, window_size = 200, is_sample = TRUE){
  #-----------------------------------------------------------------------------
  # seq_kpss: sequential KPSS test
  # 
  # Arguments:
  #   x: numeric (or vector)
  #   window_size: int
  #
  # Output:
  #   res: int (the first iteration from when the p-value is below 0.05 for 10 consequtive windows)
  #-----------------------------------------------------------------------------
  
  #dependencies
  require(tseries)
  
  n = length(x) / window_size
  stationarity = rep(0, n)
  
  #run kpss for every non-overlapping subset of length window size
  # count up the number of consequtive windows where the test is non-significant
  for(i in 1:n){
    x_sub = x[((i-1)*window_size + 1):(i*window_size)]
    stationarity[i] = as.integer(kpss.test(x_sub)$p.value > 0.05)
  }
  
  return(stationarity)
}

get_baseline = function(res, window_size = 200){
  #-----------------------------------------------------------------------------
  # get_baselineMSE: calculate the baseline MSE of a simulation
  # 
  # Arguments:
  #   res: matrix (the result of a simulation)
  #   window_size: int (the window size where we average over the results per person)
  #   is_system: logical (if it is system or person level)
  #
  # Output:
  #   MSE_baseline: matrix[numeric, nrow = nrow(res), ncol = 2 OR 3] (the mean MSE of the simulation)
  #-----------------------------------------------------------------------------
  
  #save K and theta values to the output matrix
  bv = rowMeans(res[,(ncol(res)-window_size+1):ncol(res)])
  
  return(bv)
}