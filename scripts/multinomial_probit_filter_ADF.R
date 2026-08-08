ADF_prediction <- function(a_tm1, P_tm1, parameters)
{
  return(list(a_t_tm1 = parameters$G %*% a_tm1, 
              P_t_tm1 = parameters$G %*% P_tm1 %*% t(parameters$G) + parameters$W))
  
}

SUN_correction_parameters <- function(x_t, y_t, a_t_tm1, P_t_tm1, parameters)
{
  L   = dim(x_t)[3]
  m_t = dim(x_t)[1]
  
  xi_t = a_t_tm1
  Omega_t = P_t_tm1
  
  omega_t     = diag(  sqrt(diag(Omega_t)))
  inv_omega_t = diag(1/sqrt(diag(Omega_t)))
  barOmega_t  = inv_omega_t %*% Omega_t %*% inv_omega_t
  
  barX_t   = do.call(rbind, lapply(seq(1, m_t), function(i) X_minus_l(x_t[i,,], y_t[i], L)))
  Lambda_t = do.call(adiag, lapply(seq(1, m_t), function(i) V_minus_l(y_t[i], L) %*% parameters$Sigma %*% t(V_minus_l(y_t[i], L))))
  
  XOX_Lambda = barX_t %*% Omega_t %*% t(barX_t) + Lambda_t
  
  s_t     = diag(  sqrt(diag(XOX_Lambda)))
  inv_s_t = diag(1/sqrt(diag(XOX_Lambda)))
  
  Delta_t = barOmega_t %*% omega_t %*% t(barX_t) %*% inv_s_t
  
  gamma_t = inv_s_t %*% barX_t %*% a_t_tm1
  Gamma_t = inv_s_t %*% XOX_Lambda %*% inv_s_t
  
  return(list(xi_t = xi_t, 
              Omega_t = Omega_t, omega_t = omega_t, barOmega_t = barOmega_t, 
              s_t = s_t, Delta_t = Delta_t, 
              gamma_t = gamma_t, Gamma_t = Gamma_t))
}

SUN_correction_prediction_parameters <- function(x_t, a_t_tm1, P_t_tm1, parameters)
{
  L   = dim(x_t)[3]
  m_t = dim(x_t)[1]
  
  xi_t = a_t_tm1
  Omega_t = P_t_tm1
  
  omega_t     = diag(  sqrt(diag(Omega_t)))
  inv_omega_t = diag(1/sqrt(diag(Omega_t)))
  barOmega_t  = inv_omega_t %*% Omega_t %*% inv_omega_t
  
  barX_t   = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) X_minus_l(x_t[i,,], j, L)))
  Lambda_t = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) V_minus_l(j, L) %*% parameters$Sigma %*% t(V_minus_l(j, L))))
  
  
  XOX_Lambda = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) barX_t[[i]][[j]] %*% Omega_t %*% t(barX_t[[i]][[j]]) + Lambda_t[[i]][[j]]))
  
  if(L==1){
    
    s_t     = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) t(t(diag(t(t(  sqrt(diag(XOX_Lambda[[i]][[j]])))))))))
    inv_s_t = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) t(t(diag(t(t(1/sqrt(diag(XOX_Lambda[[i]][[j]])))))))))
    
  }else{
    
    s_t     = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) diag(  sqrt(diag(XOX_Lambda[[i]][[j]])))))
    inv_s_t = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) diag(1/sqrt(diag(XOX_Lambda[[i]][[j]])))))
    
  }
  
  Delta_t = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) barOmega_t %*% omega_t %*% t(barX_t[[i]][[j]]) %*% inv_s_t[[i]][[j]]))

  gamma_t = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) inv_s_t[[i]][[j]] %*% barX_t[[i]][[j]] %*% a_t_tm1))
  Gamma_t = lapply(seq(1, m_t), function(i) lapply(seq(1, L), function(j) inv_s_t[[i]][[j]] %*% XOX_Lambda[[i]][[j]] %*% inv_s_t[[i]][[j]]))
  
  return(list(gamma_t = gamma_t, Gamma_t = Gamma_t))
}

rSUN <- function(SUN_correction_parameters_list, n=1) 
{
  d = length(SUN_correction_parameters_list$xi_t)
  m = length(SUN_correction_parameters_list$gamma_t)
  
  if(d==1) {
    omega=sqrt(SUN_correction_parameters_list$Omega_t)
  } 
  else {
    omega=diag(sqrt(diag(SUN_correction_parameters_list$Omega_t)))
  }
  
  omDe = SUN_correction_parameters_list$omega_t%*%SUN_correction_parameters_list$Delta_t
  
  Gamma_t = SUN_correction_parameters_list$Gamma_t 

  invGamma = solve(Gamma_t)
  
  PsiDel = SUN_correction_parameters_list$Omega_t-omDe%*%invGamma%*%t(omDe)
  
  sampleMultNorm=t(mvtnorm::rmvnorm(n = n, mean = rep(0,d), sigma = (PsiDel+t(PsiDel))/2))
  
  if(m==1){ 
    sampleTruncNorm = TruncatedNormal::rtmvnorm(n = n, 
                                                mu = rep(0,m), sigma = (Gamma_t+t(Gamma_t))/2, 
                                                lb = -SUN_correction_parameters_list$gamma_t, ub = rep(Inf,m))
  } else{ 
    sampleTruncNorm = t(TruncatedNormal::rtmvnorm(n = n, 
                                                  mu = rep(0,m), sigma = (Gamma_t+t(Gamma_t))/2, 
                                                  lb = -SUN_correction_parameters_list$gamma_t, ub = rep(Inf,m)))
  }
  
  rep_xi = matrix(nrow = d, ncol = n)
  rep_xi[,1:n] = SUN_correction_parameters_list$xi_t
  
  sampleSUN = rep_xi + sampleMultNorm + omDe%*%invGamma%*%sampleTruncNorm
  
  return(sampleSUN)
}

ADF_sampling_correction <- function(x_t, y_t, a_t_tm1, P_t_tm1, parameters, MC_samples = 5000)
{
  L = dim(x_t)[3]
  m_t = dim(x_t)[1]
  SUN_correction_parameters_list = SUN_correction_parameters(x_t, y_t, a_t_tm1, P_t_tm1, parameters)
  
  log_likelihood_increment = log(Phi(SUN_correction_parameters_list$gamma_t, 
                                     SUN_correction_parameters_list$Gamma_t))
  
  SUN_correction_prediction_parameters_list = SUN_correction_prediction_parameters(x_t, a_t_tm1, P_t_tm1, parameters)
  # prediction_prob = sapply(seq(1, m_t), function(i) sapply(seq(1, L), 
  #                                                          function(j)  Phi(SUN_correction_prediction_parameters_list$gamma_t[[i]][[j]], 
  #                                                                           SUN_correction_prediction_parameters_list$Gamma_t[[i]][[j]])))
  # 

  SUN_sample = rSUN(SUN_correction_parameters_list, MC_samples) 

  a_t = t(t(apply(SUN_sample, 1, mean)))
  P_t = cov(t(SUN_sample))
  
  return(list(log_likelihood_increment = log_likelihood_increment,
              # prediction_prob = prediction_prob,
              a_t = a_t,
              P_t = P_t,
              sample = SUN_sample))
}


ADF_sampling_filter <- function(data, parameters, MC_samples = 1000)
{
  start_time <- Sys.time()
  L = dim(parameters$Sigma)[1]
  Time = length(data$Y)
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  sampling_history = array(NA, dim = c(Time, dim(parameters$a_0)[1], MC_samples))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  log_likelihood = 0
  # prediction_prob_history = list() 
  
  for(t in 1:Time)
  {
    ADF_prediction_parameters = ADF_prediction(a_history[t,], P_history[t,,], parameters)
    a_t_tm1 = ADF_prediction_parameters$a_t_tm1
    P_t_tm1 = ADF_prediction_parameters$P_t_tm1
    
    y_t = data$Y[[t]]
    x_t = data$X[[t]]
    ADF_correction_output = ADF_sampling_correction(x_t, y_t, a_t_tm1, P_t_tm1, parameters, MC_samples)
    
    a_history[t+1,]  = ADF_correction_output$a_t
    P_history[t+1,,] = ADF_correction_output$P_t
    
    sampling_history[t,,] = ADF_correction_output$sample
    
    log_likelihood = log_likelihood + ADF_correction_output$log_likelihood_increment
    
    # if(L==1){
    #   prediction_prob_history[[t]] = rbind(1-ADF_correction_output$prediction_prob, 
    #                                          ADF_correction_output$prediction_prob)
    # }else{
    #   prediction_prob_history[[t]] = rbind((1-apply(ADF_correction_output$prediction_prob, 2, sum)), 
    #                                        ADF_correction_output$prediction_prob)
    # }
    
      
  }
  
  return(list(log_likelihood = log_likelihood, 
              a_history = a_history, P_history = P_history, 
              # prediction_prob_history = prediction_prob_history,
              sampling_history = sampling_history,
              running_time = as.numeric(Sys.time() - start_time, units = "secs")))
}
