gaussian_prior <- function(data, parameters)
{
  L = dim(parameters$Sigma)[1]
  Time = length(data$Y)
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  
  for(t in 1:Time)
  {
    #print(t)
    ADF_prediction_parameters = ADF_prediction(a_history[t,], P_history[t,,], parameters)
    a_t_tm1 = ADF_prediction_parameters$a_t_tm1
    P_t_tm1 = ADF_prediction_parameters$P_t_tm1
    
    
    a_history[t+1,]  = a_t_tm1
    P_history[t+1,,] = P_t_tm1
    
  }
  
  return(list(a_history = a_history, P_history = P_history))
}

multinomial_probit_EM_invSigma <- function(data, filter, initial_parameters, EM_iterations = 100, type = "std",
                                  log = list(FALSE, c(NA, NA)))
{
  Time = length(data$Y)
  sum_m_t = sum(sapply(data$Y, length)) #Reduce("+", lapply(seq(1, Time), function(t) length(data$Y[[t]])))
  
  parameters_current = initial_parameters
  
  log_likelihood_history = array(NA, dim = c(EM_iterations+1))

  invSigma_history = array(NA, dim = c(EM_iterations+1, dim(parameters_current$Sigma)[1], dim(parameters_current$Sigma)[1]))
  
  for(iteration in 1:EM_iterations)
  {
    print("##########################################################")
    print(paste("Iteration nr ", iteration))
    print("Matrix invSigma")
    print(parameters_current$Sigma)
    invSigma_history[iteration,,] = solve(parameters_current$Sigma)
    
    filter_output = filter(data, parameters_current)
    smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
    log_likelihood_history[iteration] = filter_output$log_likelihood
    print(filter_output$log_likelihood)

    eps_history = epsilon_smoother(data, parameters_current, filter_output, smoother_output, type)
      
    eps_SUN_sample = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) rSUN(eps_history[[t]][[i]], 1000)))
      
    eps_second_moment = lapply(seq(1, Time),
                                 function(t) lapply(seq(1, length(data$Y[[t]])),
                                                    function(i) cov(t(eps_SUN_sample[[t]][[i]])) + (t(t(apply(eps_SUN_sample[[t]][[i]], 1, mean))) %*% t(apply(eps_SUN_sample[[t]][[i]], 1, mean)))))
      
    new_Sigma = Reduce("+", lapply(seq(1, Time), function(t) Reduce("+", eps_second_moment[[t]])))/sum_m_t
      
    parameters_current$Sigma = (new_Sigma + t(new_Sigma))/2
    
    if(log[[1]])
    {
      log_likelihood_history_subset = log_likelihood_history[1:iteration]
      Sigma_history_subset = Sigma_history[1:iteration, , ]
      
      EM_subset_output = list( log_likelihood_history = log_likelihood_history_subset, 
                               Sigma_history = Sigma_history_subset)
      
      save(EM_subset_output, file = paste(paste(log[[2]][1], sep = "", log[[2]][2]), sep = "", "_subset_EM.Rdata"))
      
    }
    
  }
  
  print("##########################################################")
  print("END")
  print("Matrix Sigma")
  print(parameters_current$Sigma)
  Sigma_history[iteration,,] = parameters_current$Sigma
  
  filter_output = filter(data, parameters_current)
  smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
  log_likelihood_history[iteration+1] = filter_output$log_likelihood
  print(filter_output$log_likelihood)
  
  return(list(log_likelihood_history = log_likelihood_history, 
              G_history = G_history, 
              W_history = W_history,
              Sigma_history = Sigma_history))
}
