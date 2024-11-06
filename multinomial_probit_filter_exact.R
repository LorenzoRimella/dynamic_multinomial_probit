library(mvtnorm)
library(TruncatedNormal)
library(magic)

library(MomTrunc)
library(tlrmvnmvt)

exact_prediction <- function(SUN_parameters_t, parameters)
{
  xi_next = parameters$G %*% SUN_parameters_t$xi_t
    
  Omega_next = parameters$G %*% SUN_parameters_t$Omega_t %*% t(parameters$G) + parameters$W
  
  omega_next     = diag(sqrt(diag(Omega_next)))
  inv_omega_next = diag(1/sqrt(diag(Omega_next)))
  barOmega_next  = inv_omega_next %*% Omega_next %*% inv_omega_next
    
  Delta_next = inv_omega_next %*% parameters$G %*% SUN_parameters_t$omega_t %*% SUN_parameters_t$Delta_t
  
  gamma_next = SUN_parameters_t$gamma_t
    
  Gamma_bext = SUN_parameters_t$Gamma_t
  
  return(list(xi_t_tm1 = xi_next, 
              Omega_t_tm1 = Omega_next, omega_t_tm1 = omega_next, barOmega_t_tm1 = barOmega_next, 
              Delta_t_tm1 = Delta_next, 
              gamma_t_tm1 = gamma_next, Gamma_t_tm1 = Gamma_bext))
  
}

exact_correction <- function(x_t, y_t, SUN_prediction_parameters, parameters)
{
  L   = dim(x_t)[3]
  m_t = dim(x_t)[1]
  
  xi_t    = SUN_prediction_parameters$xi_t_tm1
  Omega_t = SUN_prediction_parameters$Omega_t_tm1
  
  omega_t = SUN_prediction_parameters$omega_t_tm1
  inv_omega_t = diag(1/sqrt(diag(Omega_t)))
  barOmega_t  = inv_omega_t %*% Omega_t %*% inv_omega_t
  
  barX_t   = do.call(rbind, lapply(seq(1, m_t), function(i) X_minus_l(x_t[i,,], y_t[i], L)))
  Lambda_t = do.call(adiag, lapply(seq(1, m_t), function(i) V_minus_l(y_t[i], L) %*% parameters$Sigma %*% t(V_minus_l(y_t[i], L))))
  
  XOX_Lambda = barX_t %*% Omega_t %*% t(barX_t) + Lambda_t
  
  s_t = diag(sqrt(diag(XOX_Lambda)))
  inv_s_t = diag(1/sqrt(diag(XOX_Lambda)))
  
  Delta_t = cbind(SUN_prediction_parameters$Delta_t_tm1, barOmega_t %*% omega_t %*% t(barX_t) %*% inv_s_t)
  
  gamma_t = rbind(SUN_prediction_parameters$gamma_t_tm1, inv_s_t %*% barX_t %*% xi_t)
  
  Gamma_t_r1 = cbind(SUN_prediction_parameters$Gamma_t_tm1, t(SUN_prediction_parameters$Delta_t_tm1) %*% omega_t %*% t(barX_t) %*% inv_s_t)
  Gamma_t_r2 = cbind(inv_s_t %*% barX_t %*% omega_t %*% SUN_prediction_parameters$Delta_t_tm1, inv_s_t %*% XOX_Lambda %*% inv_s_t)
  Gamma_t = rbind(Gamma_t_r1, Gamma_t_r2)
  
  return(list(xi_t = xi_t, 
              Omega_t = Omega_t, omega_t = omega_t, barOmega_t = barOmega_t, 
              s_t = s_t, Delta_t = Delta_t, 
              gamma_t = gamma_t, Gamma_t = Gamma_t))
}


exact_filter <- function(data, parameters, sliding_window, MC_samples = 100)
{
  Time = length(data$Y)
  L = dim(parameters$Sigma)[1]
  
  SUN_parameters_history = list()
  
  SUN_prediction_parameters = ADF_prediction(parameters$a_0, parameters$P_0, parameters)
  
  y_t = data$Y[[1]]
  x_t = data$X[[1]]
  SUN_parameters_history[[1]] = SUN_correction_parameters(x_t, y_t, SUN_prediction_parameters$a_t_tm1, SUN_prediction_parameters$P_t_tm1, parameters)
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  
  log_likelihood = 0
  
  for(t in 2:Time)
  {
    #print(t-1)
    SUN_sample = rSUN(SUN_parameters_history[[t-1]], MC_samples) 
    
    a_t = t(t(apply(SUN_sample, 1, mean)))
    P_t = cov(t(SUN_sample))
    
    a_history[t,]  = a_t
    P_history[t,,] = P_t
    
    if(t%%sliding_window==0){
      #print("new window")
      SUN_parameters_history[[t-1]] = SUN_correction_parameters(data$X[[t-1]], data$Y[[t-1]], a_t, P_t, parameters)
    }
    
    SUN_prediction_parameters = exact_prediction(SUN_parameters_history[[t-1]], parameters)
    
    y_t = data$Y[[t]]
    x_t = data$X[[t]]
    SUN_correction_parameters_output = exact_correction(x_t, y_t, SUN_prediction_parameters, parameters)
    
    SUN_parameters_history[[t]] = SUN_correction_parameters_output
    
    log_likelihood = log_likelihood + log(Phi(SUN_correction_parameters_output$gamma_t, 
                                              SUN_correction_parameters_output$Gamma_t))
    
  }
  #print(t)
  SUN_sample = rSUN(SUN_parameters_history[[Time]], MC_samples) 
  
  a_t = t(t(apply(SUN_sample, 1, mean)))
  P_t = cov(t(SUN_sample))
  
  a_history[Time+1,]  = a_t
  P_history[Time+1,,] = P_t
  
  return(list(log_likelihood = log_likelihood, a_history = a_history, P_history = P_history))
}