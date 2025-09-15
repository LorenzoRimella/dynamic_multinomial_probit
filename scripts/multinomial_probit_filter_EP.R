library(mvtnorm)
library(TruncatedNormal)
library(magic)

library(MomTrunc)
library(tlrmvnmvt)

SUN_site_i_parameters <- function(i, x_t, y_t, r_minus_i, Q_minus_i, parameters)
{
  
  L   = dim(x_t)[3]
  
  P_minus_i = solve(Q_minus_i)
  P_minus_i = (P_minus_i + t(P_minus_i))/2
  a_minus_i = P_minus_i %*% r_minus_i
    
  xi = a_minus_i
  Omega = P_minus_i
  
  omega     = diag(  sqrt(diag(Omega)))
  inv_omega = diag(1/sqrt(diag(Omega)))
  barOmega  = inv_omega %*% Omega %*% inv_omega
  
  barX   = X_minus_l(x_t[i-1,,], y_t[i-1], L)
  Lambda = V_minus_l(y_t[i-1], L) %*% parameters$Sigma %*% t(V_minus_l(y_t[i-1], L))
  XOX_Lambda = barX %*% Omega %*% t(barX) + Lambda
  
  if(L==1){
    s     = t(t(diag(t(t(  sqrt(diag(XOX_Lambda)))))))
    inv_s = t(t(diag(t(t(1/sqrt(diag(XOX_Lambda))))))) # double transpose cause there was a dim bug
  }else{
    s     = diag(  sqrt(diag(XOX_Lambda)))
    inv_s = diag(1/sqrt(diag(XOX_Lambda)))
  }
  
  Delta = barOmega %*% omega %*% t(barX) %*% inv_s
  
  gamma = inv_s %*% barX %*% a_minus_i
  Gamma = inv_s %*% XOX_Lambda %*% inv_s
  Gamma = (Gamma + t(Gamma))/2
  
  return(list(xi_t = xi, 
              Omega_t = Omega, omega_t = omega, barOmega_t = barOmega, 
              s_t = s, Delta_t = Delta, 
              gamma_t = gamma, Gamma_t = Gamma))
}

log_Psi<- function(r_Psi, Q_Psi)
{
  p = dim(Q_Psi)[1]
  
  return( (1/2)*(t(r_Psi)%*% solve(Q_Psi) %*% r_Psi) + (p/2)*log(2*pi) - (1/2)*log(abs(det(Q_Psi))))
}

# stabilize EP by convex combination alpha=0.999 (close to new update)

site_i_update <- function(i, x_t, y_t, r, Q, logZ, parameters, MC_samples=1e4)
{
  r_minus_i = t(t(apply(r[-i,], 2, sum)))
  Q_minus_i = apply(Q[-i,,], MARGIN = c(2,3), sum)
  
  Q_minus_i = (Q_minus_i+t(Q_minus_i))/2
  
  if(all(eigen(Q_minus_i, symmetric = TRUE, only.values = TRUE)$values>=0)){
    
    SUN_parameters_list_i = SUN_site_i_parameters(i, x_t, y_t, r_minus_i, Q_minus_i, parameters)
    
    SUN_sample = rSUN(SUN_parameters_list_i, MC_samples) 
    mu_i = t(t(apply(SUN_sample, 1, mean)))
    Sigma_i = cov(t(SUN_sample))
    
    Q[i,,] = solve(Sigma_i) - Q_minus_i
    r[i,]  = (Q_minus_i + Q[i,,]) %*% mu_i - r_minus_i
    
    r_EP = t(t(apply(r, 2, sum)))
    Q_EP = apply(Q, MARGIN = c(2,3), sum)
    
    logZ_i = log_Psi(r_EP, Q_EP) - log_Psi(r_minus_i, Q_minus_i) - log(Phi(SUN_parameters_list_i$gamma_t, SUN_parameters_list_i$Gamma_t))
    logZ[i-1]  = logZ_i
  }
  
  return(list(r = r, Q = Q, logZ = logZ))
}

EP_log_likelihood_increments <- function(r, Q, logZ)
{
  Q_0 = Q[1,,]
  r_0 = r[1,]
  
  r_sum = t(t(apply(r, 2, sum)))
  Q_sum = apply(Q, MARGIN = c(2,3), sum)
  
  log_likelihood_increment = log_Psi(r_sum, Q_sum) - log_Psi(r_0, Q_0) - sum(logZ)
  
  return(log_likelihood_increment)
}

EP_sampling_correction <- function(x_t, y_t, a_t_tm1, P_t_tm1, parameters, EP_max_iterations=1e3, MC_samples=1e4, tol)
{
  L = dim(x_t)[3]
  m_t = dim(x_t)[1]
  
  Q_0 = solve(P_t_tm1)
  r_0 = Q_0 %*% a_t_tm1
  
  Q = array(0, dim = c(m_t+1, length(r_0), length(r_0)))
  r = array(0, dim = c(m_t+1, length(r_0)))
  logZ = array(0, dim = c(m_t))
  
  Q[1,,] = Q_0
  r[1,]  = r_0
  
  Q_current = Q
  r_current = r
  logZ_current = logZ
  iteration = 1
  max_diff_vec = array(NA, dim = c(EP_max_iterations+1))
  max_diff= 100
  while((iteration<EP_max_iterations)&(max_diff>tol))
  {
    for(i in 2:(m_t+1))
    { 
      EP_parameters = site_i_update(i, x_t, y_t, r, Q, logZ, parameters, MC_samples)
      
      Q = EP_parameters$Q
      r = EP_parameters$r
      
      logZ = EP_parameters$logZ
    } 
    
    new_EP_loglike = EP_log_likelihood_increments(r, Q, logZ)
    old_EP_loglike = EP_log_likelihood_increments(r_current, Q_current, logZ_current)
    max_diff_vec[iteration] = abs( new_EP_loglike - old_EP_loglike)/abs(old_EP_loglike)
    
    max_diff = max_diff_vec[iteration]
    
    Q_current = Q
    r_current = r
    logZ_current = logZ
    iteration = iteration+1
  }
  
  r_sum = t(t(apply(r, 2, sum)))
  Q_sum = apply(Q, MARGIN = c(2,3), sum)

  P_t = solve(Q_sum)
  P_t = (P_t + t(P_t))/2
  a_t = P_t %*% r_sum
  
  log_likelihood_increment = log_Psi(r_sum, Q_sum) - log_Psi(r_0, Q_0) - sum(logZ)
  
  return(list(log_likelihood_increment = log_likelihood_increment,
              a_t = a_t,
              P_t = P_t,
              r_t = r,
              Q_t = Q))
}


EP_sampling_filter <- function(data, parameters, EP_max_iterations=1e3, MC_samples=1e4, tol = 10e-1)
{
  start_time <- Sys.time()
  L = dim(parameters$Sigma)[1]
  Time = length(data$Y)
  
  V_minus_pre_computed = lapply(seq(0, L), function(i) V_minus_l(i,L))
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  sampling_history = array(NA, dim = c(Time, dim(parameters$a_0)[1], MC_samples))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  log_likelihood = 0
  
  r_history = list() 
  Q_history = list() 
  
  for(t in 1:Time)
  {
    # print(t)
    EP_prediction_parameters = ADF_prediction(a_history[t,], P_history[t,,], parameters)
    a_t_tm1 = EP_prediction_parameters$a_t_tm1
    P_t_tm1 = EP_prediction_parameters$P_t_tm1
    
    y_t = data$Y[[t]]
    x_t = data$X[[t]]
    EP_correction_output = EP_sampling_correction(x_t, y_t, a_t_tm1, P_t_tm1, parameters, EP_max_iterations, MC_samples, tol)
    
    a_history[t+1,]  = EP_correction_output$a_t
    P_history[t+1,,] = EP_correction_output$P_t 
    
    sampling_history[t,,] = mvtnorm::rmvnorm(n = MC_samples, mean = EP_correction_output$a_t, sigma = EP_correction_output$P_t)
    
    log_likelihood = log_likelihood + EP_correction_output$log_likelihood_increment
    
    r_history[[t]] = EP_correction_output$r_t
    Q_history[[t]] = EP_correction_output$Q_t
    
  }
  
  return(list(log_likelihood = log_likelihood, 
              a_history = a_history, P_history = P_history, 
              r_history = r_history, Q_history = Q_history,
              sampling_history = sampling_history,
              running_time = as.numeric(Sys.time() - start_time, units = "secs")))
}

EP_predictive <- function(data, parameters, EP_output, MC_samples)
{
  L = dim(parameters$Sigma)[1]
  Time = dim(EP_output$a_history)[1]-1
  
  V_minus_pre_computed = lapply(seq(0, L), function(i) V_minus_l(i,L))
  
  single_t <- function(t)
  {
    # print(t)
    EP_prediction_parameters = ADF_prediction(EP_output$a_history[t,], EP_output$P_history[t,,], parameters)
    
    sample_beta_t_tm1 = mvtnorm::rmvnorm(n = MC_samples, mean = EP_prediction_parameters$a_t_tm1, sigma = EP_prediction_parameters$P_t_tm1)
    
    
    x_t = data$X[[t]]
    
    utility_i <- function(i) t(t(x_t[i,,]) %*% t(sample_beta_t_tm1)) + mvtnorm::rmvnorm(n = MC_samples, mean = t(t(rep(0, L))), sigma = parameters$Sigma)
    
    return(t(sapply(seq(1, dim(x_t)[1]), function(i) as.numeric(table(factor(apply(cbind(rep(0, L), utility_i(i)), 1, which.max), levels = 1:L)))/MC_samples)))
  }
  
  return(lapply(seq(1, Time), single_t))
}

