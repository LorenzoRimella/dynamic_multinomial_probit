SMC_sampling_filter <- function(data, parameters, MC_samples)
{
  
  start_time <- Sys.time()
  L = dim(parameters$Sigma)[1]
  Time = length(data$Y)
  
  V_minus_pre_computed = lapply(seq(0, L), function(i) V_minus_l(i,L))
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  sampling_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], MC_samples))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  log_likelihood = 0
  
  sampling_history[1,,] = t(mvrnorm(n = MC_samples, parameters$a_0, parameters$P_0)) 
  ESS_history = list()
  
  for(t in 1:Time)
  {
    sampling_prediction = parameters$G %*% sampling_history[t,,] + t(mvrnorm(n = MC_samples, 0*parameters$a_0, parameters$W))
    
    y_t = data$Y[[t]]
    x_t = data$X[[t]]
    
    m_t = dim(x_t)[1]
    
    # prod_likelihood <-function(s)
    # {
    #   prod(sapply(seq(1, m_t), function(i) Phi(t(t((X_minus_l(x_t[i,,], y_t[i], L) %*% sampling_prediction)[,s])),
    #                                            V_minus_l(y_t[i], L) %*% parameters$Sigma %*% t(V_minus_l(y_t[i], L)))))
    # }
    
    # start_time <- Sys.time()
    before_weights = sapply(seq(1, m_t), function(i) Phi_batch(t(X_minus_l(x_t[i,,], y_t[i], L) %*% sampling_prediction[,1:MC_samples]),
                                                 V_minus_l(y_t[i], L) %*% parameters$Sigma %*% t(V_minus_l(y_t[i], L))))
    weights = apply(before_weights, 1, prod)
    # print(as.numeric(Sys.time() - start_time))
    
    # start_time <- Sys.time()
    # weights = sapply(seq(1, MC_samples), function(s) prod_likelihood(s))
    # print(as.numeric(Sys.time() - start_time))
    
    weights_norm = weights/sum(weights)
    
    log_likelihood = log_likelihood + log(mean(weights))
    
    ancestors = sample(seq(1, MC_samples), MC_samples, replace = TRUE, prob = weights_norm)
    
    sampling_history[t+1,,] = sampling_prediction[,ancestors]
    
    a_t = t(t(apply(sampling_history[t+1,,], 1, mean)))
    P_t = cov(t(sampling_history[t+1,,]))
    
    a_history[t+1,]  = a_t
    P_history[t+1,,] = P_t
    
    ESS_history[[t]] = 1/sum(weights_norm*weights_norm)
    
  }
  
  return(list(log_likelihood = log_likelihood, 
              a_history = a_history, P_history = P_history,
              sampling_history = sampling_history[2:(Time+1),,],
              running_time = as.numeric(Sys.time() - start_time, units = "secs"),
              ESS_history = ESS_history))
}