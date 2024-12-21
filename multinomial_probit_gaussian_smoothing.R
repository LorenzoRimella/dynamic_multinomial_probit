SUN_translator_eps <- function(x_t, y_t, a, P, parameters)
{
  L   = dim(x_t)[3]
  m_t = dim(x_t)[1]
  
  xi0_t = rep(0, dim(parameters$Sigma)[1])
  Omega0_t = parameters$Sigma
  
  omega0_t = diag(sqrt(diag(Omega0_t)))
  inv_omega0_t = diag(1/sqrt(diag(Omega0_t)))
  barOmega0_t  = inv_omega0_t %*% Omega0_t %*% inv_omega0_t
  
  XPX_VSigV = lapply(seq(1, m_t), function(i) X_minus_l(x_t[i,,], y_t[i], L) %*% P[[i]] %*% t(X_minus_l(x_t[i,,], y_t[i], L)) + (-V_minus_l(y_t[i], L)) %*% Omega0_t %*% t(-V_minus_l(y_t[i], L)) )
  
  s0_t     = lapply(seq(1, m_t), function(i) diag(sqrt(diag(XPX_VSigV[[i]]))))
  inv_s0_t = lapply(seq(1, m_t), function(i) diag(1/sqrt(diag(XPX_VSigV[[i]]))))
  
  Delta0_t = lapply(seq(1, m_t), function(i) barOmega0_t %*% omega0_t %*% t(-V_minus_l(y_t[i], L)) %*% inv_s0_t[[i]])
  
  gamma0_t = lapply(seq(1, m_t), function(i) inv_s0_t[[i]] %*% X_minus_l(x_t[i,,], y_t[i], L) %*% a[[i]])
  Gamma0_t = lapply(seq(1, m_t), function(i) inv_s0_t[[i]] %*% XPX_VSigV[[i]] %*% inv_s0_t[[i]])
  
  output_SUN = lapply(seq(1, m_t), function(i) list(xi_t = xi0_t, 
                                                    Omega_t = Omega0_t, omega_t = omega0_t, barOmega_t = barOmega0_t, 
                                                    s_t = s0_t[[i]], Delta_t = Delta0_t[[i]], 
                                                    gamma_t = gamma0_t[[i]], Gamma_t = Gamma0_t[[i]]))
  
  return(output_SUN)
}

Gaussian_smoother <- function(data, parameters, Gaussian_filter)
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
    P_tp1_t = (P_tp1_t +t(P_tp1_t))/2
    J_t     = P_t %*% t(parameters$G) %*% solve(P_tp1_t)
    
    a_smoothing_history[Time-t+1,]  = a_t + J_t %*% (a_smoothing_history[Time-t+2,] - a_tp1_t)
    P_smoothing_history[Time-t+1,,] = J_t %*% P_smoothing_history[Time-t+2,,] %*% t(J_t) + P_t - J_t %*% parameters$G %*% P_t
    P_smoothing_history[Time-t+1,,] = (P_smoothing_history[Time-t+1,,] + t(P_smoothing_history[Time-t+1,,]))/2
    J_history[Time-t+1,,] = J_t
  }
  
  return(list(a_history = a_smoothing_history, P_history = P_smoothing_history, J_history = J_history))
}

epsilon_smoother <- function(data, parameters, EP_filter_output, EP_smoother_output)
{
  Time = dim(EP_filter_output$a_history)[1]-1
  
  a_minus_i_smoothing_history = list()
  P_minus_i_smoothing_history = list()
  
  eps_smoothing_history = list()
  
  for(t in 1:(Time))
  {
    # print(t)
    y_t = data$Y[[t]]
    x_t = data$X[[t]]
    
    m_t = dim(x_t)[1]
    
    Q_t_T = solve(EP_smoother_output$P_history[t+1,,])
    r_t_T = Q_t_T %*% EP_smoother_output$a_history[t+1,]

    Q_t_minus_i_T = lapply(seq(1, m_t), function(i) Q_t_T - EP_filter_output$Q[[t]][i+1,,])
    r_t_minus_i_T = lapply(seq(1, m_t), function(i) r_t_T - t(t(EP_filter_output$r[[t]][i+1,])))

    P_t_minus_i_T = lapply(seq(1, m_t), function(i) solve(Q_t_minus_i_T[[i]]))
    a_t_minus_i_T = lapply(seq(1, m_t), function(i) P_t_minus_i_T[[i]] %*% r_t_minus_i_T[[i]])
    
    # P_t_minus_i_T = lapply(seq(1, m_t), function(i) EP_smoother_output$P_history[t+1,,])
    # a_t_minus_i_T = lapply(seq(1, m_t), function(i) EP_smoother_output$a_history[t+1,])
    
    eps_smoothing_t = SUN_translator_eps(x_t, y_t, a_t_minus_i_T, P_t_minus_i_T, parameters)
    
    eps_smoothing_history[[t]] = eps_smoothing_t
  }
  
  return(eps_smoothing_history)
  
}


# a_minus_i <- function(i, r, Q)
# {
#   r_minus_i = t(t(apply(r[-i,], 2, sum)))
#   Q_minus_i = apply(Q[-i,,], MARGIN = c(2,3), sum)
#   
#   P_minus_i = solve(Q_minus_i)
#   P_minus_i = (P_minus_i + t(P_minus_i))/2
#   a_minus_i = P_minus_i %*% r_minus_i
#   
#   return(a_minus_i)
# }
# 
# P_minus_i <- function(i, r, Q)
# {
#   r_minus_i = t(t(apply(r[-i,], 2, sum)))
#   Q_minus_i = apply(Q[-i,,], MARGIN = c(2,3), sum)
#   
#   P_minus_i = solve(Q_minus_i)
#   P_minus_i = (P_minus_i + t(P_minus_i))/2
#   
#   return(P_minus_i)
# }

# epsilon_smoother <- function(data, parameters, EP_filter_output, EP_smoother_output)
# {
#   Time = dim(EP_filter_output$a_history)[1]-1
#   
#   a_minus_i_smoothing_history = list()
#   P_minus_i_smoothing_history = list()
#   
#   eps_smoothing_history = list()
#   
#   for(t in 1:(Time-1))
#   {
#     # print(t)
#     y_t = data$Y[[t]]
#     x_t = data$X[[t]]
#     
#     m_t = dim(x_t)[1]
#     
#     a_minus_i_t = lapply(seq(1, m_t), function(i) a_minus_i(i, EP_filter_output$r_history[[t]], EP_filter_output$Q_history[[t]]))
#     P_minus_i_t = lapply(seq(1, m_t), function(i) P_minus_i(i, EP_filter_output$r_history[[t]], EP_filter_output$Q_history[[t]]))
#     
#     a_tp1_minus_i_t = lapply(seq(1, m_t), function(i) parameters$G %*% a_minus_i_t[[i]])
#     P_tp1_minus_i_t = lapply(seq(1, m_t), function(i) parameters$G %*% P_minus_i_t[[i]] %*% t(parameters$G) + parameters$W)
#     
#     J_minus_i_t     = lapply(seq(1, m_t), function(i) P_minus_i_t[[i]] %*% t(parameters$G) %*% solve(P_tp1_minus_i_t[[i]]))
#     
#     a_t_minus_i_T = lapply(seq(1, m_t), function(i) a_minus_i_t[[i]] + J_minus_i_t[[i]] %*% (EP_smoother_output$a_history[t+2,] - a_tp1_minus_i_t[[i]]))
#     P_t_minus_i_T = lapply(seq(1, m_t), function(i) J_minus_i_t[[i]] %*% EP_smoother_output$P_history[t+2,,] %*% t(J_minus_i_t[[i]]) + P_minus_i_t[[i]] - J_minus_i_t[[i]] %*% parameters$G %*% P_minus_i_t[[i]])
#     P_t_minus_i_T = lapply(seq(1, m_t), function(i) (P_t_minus_i_T[[i]]+t(P_t_minus_i_T[[i]]))/2)
#     
#     # eps_smoothing_t = SUN_translator_eps(x_t, y_t, a_t_minus_i_T, P_t_minus_i_T, parameters)
#     eps_smoothing_t = SUN_translator_eps(x_t, y_t, a_minus_i_t, P_minus_i_t, parameters)
#     
#     eps_smoothing_history[[t]] = eps_smoothing_t
#   }
#   
#   t = Time
#   
#   y_t = data$Y[[t]]
#   x_t = data$X[[t]]
#   
#   m_t = dim(x_t)[1]
#   
#   a_t_minus_i_T = lapply(seq(1, m_t), function(i) a_minus_i(i, EP_filter_output$r_history[[t]], EP_filter_output$Q_history[[t]]))
#   P_t_minus_i_T = lapply(seq(1, m_t), function(i) P_minus_i(i, EP_filter_output$r_history[[t]], EP_filter_output$Q_history[[t]]))
#   P_t_minus_i_T = lapply(seq(1, m_t), function(i) (P_t_minus_i_T[[i]]+t(P_t_minus_i_T[[i]]))/2)
#   
#   # eps_smoothing_t = SUN_translator_eps(x_t, y_t, a_t_minus_i_T, P_t_minus_i_T, parameters)
#   eps_smoothing_t = SUN_translator_eps(x_t, y_t, a_minus_i_t, P_minus_i_t, parameters)
#   
#   eps_smoothing_history[[t]] = eps_smoothing_t
#   
#   return(eps_smoothing_history)
#   
# }

