Gaussian_smoother <- function(parameters, Gaussian_filter)
{
  Time = dim(Gaussian_filter$a_history)[1]-1
  
  a_smoothing_history = array(NA, dim = c(Time+1, dim(Gaussian_filter$a_history)[2]))
  P_smoothing_history = array(NA, dim = c(Time+1, dim(Gaussian_filter$a_history)[2], dim(Gaussian_filter$a_history)[2]))
  J_history = array(NA, dim = c(Time, dim(Gaussian_filter$a_history)[2], dim(Gaussian_filter$a_history)[2]))
  
  a_smoothing_history[Time+1,]  = Gaussian_filter$a_history[Time+1,]
  P_smoothing_history[Time+1,,] = Gaussian_filter$P_history[Time+1,,]
  
  for(t in 1:Time)
  {
    a_t = Gaussian_filter$a_history[Time-t+1,]
    P_t = Gaussian_filter$P_history[Time-t+1,,]
    
    a_tp1_t = parameters$G %*% a_t
    P_tp1_t = parameters$G %*% P_t %*% t(parameters$G) + parameters$W
    J_t     = P_t %*% t(parameters$G) %*% solve(P_tp1_t)
    
    a_smoothing_history[Time-t+1,]  = a_t + J_t %*% (a_smoothing_history[Time-t+2,] - a_tp1_t)
    P_smoothing_history[Time-t+1,,] = J_t %*% P_smoothing_history[Time-t+2,,] %*% t(J_t) + P_t - J_t %*% parameters$G %*% P_t
    J_history[Time-t+1,,] = J_t
  }
  
  return(list(a_history = a_smoothing_history, P_history = P_smoothing_history, J_history = J_history))
}
