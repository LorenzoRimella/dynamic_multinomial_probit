SUN_translator <- function(x_t, y_t, a_t_T, P_t_T, parameters)
{
  L   = dim(x_t)[3]
  m_t = dim(x_t)[1]
  
  xi0_t = rep(0, dim(parameters$Sigma)[1])
  Omega0_t = parameters$Sigma
  
  omega0_t = diag(sqrt(diag(Omega0_t)))
  inv_omega0_t = diag(1/sqrt(diag(Omega0_t)))
  barOmega0_t  = inv_omega0_t %*% Omega0_t %*% inv_omega0_t
  
  XPX_VSigV = lapply(seq(1, m_t), function(i) X_minus_l(x_t[i,,], y_t[i], L) %*% P_t_T %*% t(X_minus_l(x_t[i,,], y_t[i], L)) + (-V_minus_l(y_t[i], L)) %*% Omega0_t %*% t(-V_minus_l(y_t[i], L)) )
  
  s0_t     = lapply(seq(1, m_t), function(i) diag(sqrt(diag(XPX_VSigV[[i]]))))
  inv_s0_t = lapply(seq(1, m_t), function(i) diag(1/sqrt(diag(XPX_VSigV[[i]]))))
  
  Delta0_t = lapply(seq(1, m_t), function(i) barOmega0_t %*% omega0_t %*% t(-V_minus_l(y_t[i], L)) %*% inv_s0_t[[i]])
  
  #XPX_VDeltV = lapply(seq(1, m_t), function(i) X_minus_l(x_t[i,,], y_t[i], L) %*% P_t_T %*% t(X_minus_l(x_t[i,,], y_t[i], L)) + (-V_minus_l(y_t[i], L)) %*% Delta0_t[[i]] %*% t(-V_minus_l(y_t[i], L)) )
  gamma0_t = lapply(seq(1, m_t), function(i) inv_s0_t[[i]] %*% X_minus_l(x_t[i,,], y_t[i], L) %*% a_t_T)
  Gamma0_t = lapply(seq(1, m_t), function(i) inv_s0_t[[i]] %*% XPX_VSigV[[i]] %*% inv_s0_t[[i]])
  
  output_SUN = lapply(seq(1, m_t), function(i) list(xi_t = xi0_t, 
                                                    Omega_t = Omega0_t, omega_t = omega0_t, barOmega_t = barOmega0_t, 
                                                    s_t = s0_t[[i]], Delta_t = Delta0_t[[i]], 
                                                    gamma_t = gamma0_t[[i]], Gamma_t = Gamma0_t[[i]]))
  
  return(output_SUN)
}

Gaussian_smoother <- function(data, parameters, Gaussian_filter, type = "smoother")
{
  Time = dim(Gaussian_filter$a_history)[1]-1
  
  a_smoothing_history = array(NA, dim = c(Time+1, dim(Gaussian_filter$a_history)[2]))
  P_smoothing_history = array(NA, dim = c(Time+1, dim(Gaussian_filter$a_history)[2], dim(Gaussian_filter$a_history)[2]))
  J_history = array(NA, dim = c(Time, dim(Gaussian_filter$a_history)[2], dim(Gaussian_filter$a_history)[2]))
  
  a_smoothing_history[Time+1,]  = Gaussian_filter$a_history[Time+1,]
  P_smoothing_history[Time+1,,] = Gaussian_filter$P_history[Time+1,,]
  
  eps_smoothing_history = list()
  
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
    
    # Epsilons
    y_t = data$Y[[Time-t+1]]
    x_t = data$X[[Time-t+1]]
    
    if(type=="smoother"){
      eps_smoothing_t = SUN_translator(x_t, y_t, a_smoothing_history[Time-t+2,], P_smoothing_history[Time-t+2,,], parameters)
    }
    else{
      eps_smoothing_t = SUN_translator(x_t, y_t, a_tp1_t, P_tp1_t, parameters)
    }
    
    eps_smoothing_history[[Time-t+1]] = eps_smoothing_t
  }
  
  return(list(a_history = a_smoothing_history, P_history = P_smoothing_history, J_history = J_history, eps_history = eps_smoothing_history))
}
