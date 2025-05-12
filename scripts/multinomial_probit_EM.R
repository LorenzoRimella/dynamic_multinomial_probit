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

multinomial_probit_EM <- function(data, filter, initial_parameters, EM_iterations = 100, 
                                  learn_G = TRUE, learn_W = TRUE, learn_Sigma = TRUE, type = "std",
                                  log = list(FALSE, c(NA, NA)))
{
  Time = length(data$Y)
  sum_m_t = sum(sapply(data$Y, length)) #Reduce("+", lapply(seq(1, Time), function(t) length(data$Y[[t]])))
  
  parameters_current = initial_parameters
  
  log_likelihood_history = array(NA, dim = c(EM_iterations+1))
  
  G_history     = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$G)[1],     dim(parameters_current$G)[1]))
  W_history     = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$W)[1],     dim(parameters_current$W)[1]))
  Sigma_history = array(NA, dim = c(EM_iterations+1, dim(parameters_current$Sigma)[1], dim(parameters_current$Sigma)[1]))
  
  for(iteration in 1:EM_iterations)
  {
    print(paste("Iteration nr ", iteration))
    filter_output = filter(data, parameters_current)
    smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
    log_likelihood_history[iteration] = filter_output$log_likelihood
    print(filter_output$log_likelihood)
    
    if(learn_G)
    {
      print("Matrix G")
      print(parameters_current$G)
    }
    G_history[iteration,,] = parameters_current$G
    
    if(learn_W)
    {
      print("Matrix W")
      print(parameters_current$W)
    }
    W_history[iteration,,] = parameters_current$W
    
    if(learn_Sigma)
    {
      print("Matrix Sigma")
      print(parameters_current$Sigma)
    }
    Sigma_history[iteration,,] = parameters_current$Sigma
    
    if(learn_Sigma)
    {
      eps_history = epsilon_smoother(data, parameters_current, filter_output, smoother_output, type)
      
      eps_SUN_sample = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) rSUN(eps_history[[t]][[i]], 1000)))
      
      eps_second_moment = lapply(seq(1, Time),
                                 function(t) lapply(seq(1, length(data$Y[[t]])),
                                                    function(i) cov(t(eps_SUN_sample[[t]][[i]])) + (t(t(apply(eps_SUN_sample[[t]][[i]], 1, mean))) %*% t(apply(eps_SUN_sample[[t]][[i]], 1, mean)))))
      
      new_Sigma = Reduce("+", lapply(seq(1, Time), function(t) Reduce("+", eps_second_moment[[t]])))/sum_m_t
      
      parameters_current$Sigma = (new_Sigma + t(new_Sigma))/2
    }
    
    if(learn_G)
    {
      a_0_Tm1_T = smoother_output$a_history[1:Time,]
      a_1_T_T   = smoother_output$a_history[2:(Time+1),]
      
      P_0_Tm1_T = smoother_output$P_history[1:Time,,]
      P_1_T_T   = smoother_output$P_history[2:(Time+1),,]
      
      J_0_T_T   = smoother_output$J_history[1:Time,,]
      
      sum_barGP_a_t_tm1_T = Reduce("+", lapply(seq(1, Time), 
                                               function(t) J_0_T_T[t,,] %*% P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,])))
      sum_P_a_tm1_tm1_T   = Reduce("+", lapply(seq(1, Time), 
                                               function(t) P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])))
      
      new_G = sum_barGP_a_t_tm1_T %*% solve(sum_P_a_tm1_tm1_T)
      
    }else{
      new_G = parameters_current$G
    }
    
    if(learn_W)
    {
      new_W = (1/Time)*Reduce("+", lapply(seq(1, Time),
                                          function(t) P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_1_T_T[t,])
                                          - J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) - a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G)
                                          - t(J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G))
                                          + new_G %*% (P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])) %*% t(new_G)))
      
    }else{
      new_W = parameters_current$W
    }
    
    parameters_current$G = new_G
    parameters_current$W = new_W
    
    if(log[[1]])
    {
      log_likelihood_history_subset = log_likelihood_history[1:iteration]
      G_history_subset     = G_history[1:iteration, , ]
      W_history_subset     = W_history[1:iteration, , ]
      Sigma_history_subset = Sigma_history[1:iteration, , ]
      
      EM_subset_output = list( log_likelihood_history = log_likelihood_history_subset, 
                               G_history = G_history_subset, 
                               W_history = W_history_subset,
                               Sigma_history = Sigma_history_subset)
      
      save(EM_subset_output, file = paste(paste(log[[2]][1], sep = "", log[[2]][2]), sep = "", "_subset_EM.Rdata"))

    }
    
  }
  
  print(paste("END"))
  filter_output = filter(data, parameters_current)
  smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
  log_likelihood_history[iteration+1] = filter_output$log_likelihood
  print(filter_output$log_likelihood)
  
  if(learn_G)
  {
    print("Matrix G")
    print(parameters_current$G)
    G_history[iteration+1,,] = parameters_current$G
  }
  
  if(learn_W)
  {
    print("Matrix W")
    print(parameters_current$W)
    W_history[iteration+1,,] = parameters_current$W
  }
  
  if(learn_Sigma)
  {
    print("Matrix Sigma")
    print(parameters_current$Sigma)
    Sigma_history[iteration+1,,] = parameters_current$Sigma
  }
  
  return(list(log_likelihood_history = log_likelihood_history, 
              G_history = G_history, 
              W_history = W_history,
              Sigma_history = Sigma_history))
}

multinomial_probit_EM_trace <- function(data, filter, initial_parameters, EM_iterations = 100, 
                                  learn_G = TRUE, learn_W = TRUE, learn_Sigma = TRUE, type = "std",
                                  normalization = 1, log = c(FALSE, NA))
{
  Time = length(data$Y)
  sum_m_t = sum(sapply(data$Y, length)) #Reduce("+", lapply(seq(1, Time), function(t) length(data$Y[[t]])))
  L = dim(initial_parameters$Sigma)[1]
  
  parameters_current = initial_parameters
  trace_Sigma = sum(diag(parameters_current$Sigma))
  parameters_current$P_0   = normalization*parameters_current$P_0/trace_Sigma
  parameters_current$W     = normalization*parameters_current$W/trace_Sigma
  parameters_current$Sigma = normalization*parameters_current$Sigma/trace_Sigma
  
  log_likelihood_history = array(NA, dim = c(EM_iterations+1))
  
  P_0_history   = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$P_0)[1],   dim(parameters_current$P_0)[1]))
  G_history     = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$G)[1],     dim(parameters_current$G)[1]))
  W_history     = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$W)[1],     dim(parameters_current$W)[1]))
  Sigma_history = array(NA, dim = c(EM_iterations+1, dim(parameters_current$Sigma)[1], dim(parameters_current$Sigma)[1]))
  
  for(iteration in 1:EM_iterations)
  {
    print(paste("Iteration nr ", iteration))
    filter_output = filter(data, parameters_current)
    smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
    log_likelihood_history[iteration] = filter_output$log_likelihood
    print(filter_output$log_likelihood)
    
    if(learn_G)
    {
      print("Matrix G")
      print(parameters_current$G)
    }
    G_history[iteration,,] = parameters_current$G
    
    
    if(learn_W)
    {
      print("Matrix W")
      print(parameters_current$W)
    }
    W_history[iteration,,] = parameters_current$W
    
    if(learn_Sigma)
    {
      print("Matrix Sigma")
      print(parameters_current$Sigma)
    }
    Sigma_history[iteration,,] = parameters_current$Sigma
    
    P_0_history[iteration,,] = parameters_current$P_0
    
    if(learn_Sigma)
    {
      eps_history = epsilon_smoother(data, parameters_current, filter_output, smoother_output, type)
      
      eps_SUN_sample = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) rSUN(eps_history[[t]][[i]], 1000)))
      
      eps_second_moment = lapply(seq(1, Time),
                                 function(t) lapply(seq(1, length(data$Y[[t]])),
                                                    function(i) cov(t(eps_SUN_sample[[t]][[i]])) + (t(t(apply(eps_SUN_sample[[t]][[i]], 1, mean))) %*% t(apply(eps_SUN_sample[[t]][[i]], 1, mean)))))
      
      new_Sigma = Reduce("+", lapply(seq(1, Time), function(t) Reduce("+", eps_second_moment[[t]])))/sum_m_t
      
      parameters_current$Sigma = (new_Sigma + t(new_Sigma))/2
    }
    
    if(learn_G)
    {
      a_0_Tm1_T = smoother_output$a_history[1:Time,]
      a_1_T_T   = smoother_output$a_history[2:(Time+1),]
      
      P_0_Tm1_T = smoother_output$P_history[1:Time,,]
      P_1_T_T   = smoother_output$P_history[2:(Time+1),,]
      
      J_0_T_T   = smoother_output$J_history[1:Time,,]
      
      sum_barGP_a_t_tm1_T = Reduce("+", lapply(seq(1, Time), 
                                               function(t) J_0_T_T[t,,] %*% P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,])))
      sum_P_a_tm1_tm1_T   = Reduce("+", lapply(seq(1, Time), 
                                               function(t) P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])))
      
      new_G = sum_barGP_a_t_tm1_T %*% solve(sum_P_a_tm1_tm1_T)
      
    }else{
      new_G = parameters_current$G
    }
    
    if(learn_W)
    {
      new_W = (1/Time)*Reduce("+", lapply(seq(1, Time),
                                          function(t) P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_1_T_T[t,])
                                          - J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) - a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G)
                                          - t(J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G))
                                          + new_G %*% (P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])) %*% t(new_G)))
      
    }else{
      new_W = parameters_current$W
    }
    
    parameters_current$G = new_G
    parameters_current$W = new_W
    
    trace_Sigma = sum(diag(parameters_current$Sigma))
    parameters_current$P_0   = normalization*parameters_current$P_0/trace_Sigma
    parameters_current$W     = normalization*parameters_current$W/trace_Sigma
    parameters_current$Sigma = normalization*parameters_current$Sigma/trace_Sigma
    
    if(log[[1]])
    {
      log_likelihood_history_subset = log_likelihood_history[1:iteration]
      P_0_history_subset   = P_0_history[1:iteration, , ]
      G_history_subset     = G_history[1:iteration, , ]
      W_history_subset     = W_history[1:iteration, , ]
      Sigma_history_subset = Sigma_history[1:iteration, , ]
      
      EM_subset_output = list( log_likelihood_history = log_likelihood_history_subset, 
                               P_0_history = P_0_history_subset,
                               G_history = G_history_subset, 
                               W_history = W_history_subset,
                               Sigma_history = Sigma_history_subset)
      
      save(EM_subset_output, file = paste(paste(log[[2]][1], sep = "", log[[2]][2]), sep = "", "_subset_EM.Rdata"))
      
    }
    
  }
  
  print(paste("END"))
  filter_output = filter(data, parameters_current)
  smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
  log_likelihood_history[iteration+1] = filter_output$log_likelihood
  print(filter_output$log_likelihood)
  
  if(learn_G)
  {
    print("Matrix G")
    print(parameters_current$G)
    G_history[iteration+1,,] = parameters_current$G
  }
  
  if(learn_W)
  {
    print("Matrix W")
    print(parameters_current$W)
    W_history[iteration+1,,] = parameters_current$W
  }
  
  if(learn_Sigma)
  {
    print("Matrix Sigma")
    print(parameters_current$Sigma)
    Sigma_history[iteration+1,,] = parameters_current$Sigma
  }
  
  return(list(log_likelihood_history = log_likelihood_history, 
              P_0_history = P_0_history,
              G_history = G_history, 
              W_history = W_history,
              Sigma_history = Sigma_history))
}

# Function that returns A as a function of lambda
A_func <- function(lambda,sum_m_t,V) {
  L <- ncol(V)
  (1/sum_m_t) * (V - 2 * lambda * diag(1,L,L))
}

# Function to compute trace of inverse of A(c)
trace_inv_A <- function(lambda,sum_m_t,V) {
  A <- A_func(lambda,sum_m_t,V)
  
  if (det(A) <= 1e-10) {
    return(Inf)  # Avoid singular or near-singular matrices
  }
  
  A_inv <- solve(A)
  sum(diag(A_inv))
}

# Root-finding function: solve trace_inv_A(c) = L
root_function <- function(lambda,sum_m_t,V,target_trace) {
  trace_inv_A(lambda,sum_m_t,V) - target_trace
}

multinomial_probit_EM_constrained <- function(data, filter, initial_parameters, EM_iterations = 100, 
                                  learn_G = TRUE, learn_W = TRUE, learn_Sigma = TRUE, type = "std",
                                  log = list(FALSE, c(NA, NA)),
                                  target_trace)
{
  Time = length(data$Y)
  sum_m_t = sum(sapply(data$Y, length)) #Reduce("+", lapply(seq(1, Time), function(t) length(data$Y[[t]])))
  
  parameters_current = initial_parameters
  
  log_likelihood_history = array(NA, dim = c(EM_iterations+1))
  
  G_history     = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$G)[1],     dim(parameters_current$G)[1]))
  W_history     = array(NA, dim = c(EM_iterations+1,     dim(parameters_current$W)[1],     dim(parameters_current$W)[1]))
  Sigma_history = array(NA, dim = c(EM_iterations+1, dim(parameters_current$Sigma)[1], dim(parameters_current$Sigma)[1]))
  
  for(iteration in 1:EM_iterations)
  {
    print(paste("Iteration nr ", iteration))
    filter_output = filter(data, parameters_current)
    smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
    log_likelihood_history[iteration] = filter_output$log_likelihood
    print(filter_output$log_likelihood)
    
    if(learn_G)
    {
      print("Matrix G")
      print(parameters_current$G)
    }
    G_history[iteration,,] = parameters_current$G
    
    if(learn_W)
    {
      print("Matrix W")
      print(parameters_current$W)
    }
    W_history[iteration,,] = parameters_current$W
    
    if(learn_Sigma)
    {
      print("Matrix Sigma")
      print(parameters_current$Sigma)
    }
    Sigma_history[iteration,,] = parameters_current$Sigma
    
    if(learn_Sigma)
    {
      eps_history = epsilon_smoother(data, parameters_current, filter_output, smoother_output, type)
      
      eps_SUN_sample = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) rSUN(eps_history[[t]][[i]], 10000)))
      
      eps_second_moment = lapply(seq(1, Time),
                                 function(t) lapply(seq(1, length(data$Y[[t]])),
                                                    function(i) cov(t(eps_SUN_sample[[t]][[i]])) + (t(t(apply(eps_SUN_sample[[t]][[i]], 1, mean))) %*% t(apply(eps_SUN_sample[[t]][[i]], 1, mean)))))
      
      V = Reduce("+", lapply(seq(1, Time), function(t) Reduce("+", eps_second_moment[[t]])))
      
      
      result <- uniroot(root_function, lower = -10, upper = 10, 
                        sum_m_t = sum_m_t, V = V, target_trace = target_trace)
      
      lambda_opt <- result$root
      
      print(lambda_opt)
      
      new_Sigma = (1/sum_m_t) * (V - 2 * lambda_opt * diag(1,L,L))
      
      parameters_current$Sigma = (new_Sigma + t(new_Sigma))/2
    }
    
    if(learn_G)
    {
      a_0_Tm1_T = smoother_output$a_history[1:Time,]
      a_1_T_T   = smoother_output$a_history[2:(Time+1),]
      
      P_0_Tm1_T = smoother_output$P_history[1:Time,,]
      P_1_T_T   = smoother_output$P_history[2:(Time+1),,]
      
      J_0_T_T   = smoother_output$J_history[1:Time,,]
      
      sum_barGP_a_t_tm1_T = Reduce("+", lapply(seq(1, Time), 
                                               function(t) J_0_T_T[t,,] %*% P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,])))
      sum_P_a_tm1_tm1_T   = Reduce("+", lapply(seq(1, Time), 
                                               function(t) P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])))
      
      new_G = sum_barGP_a_t_tm1_T %*% solve(sum_P_a_tm1_tm1_T)
      
    }else{
      new_G = parameters_current$G
    }
    
    if(learn_W)
    {
      new_W = (1/Time)*Reduce("+", lapply(seq(1, Time),
                                          function(t) P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_1_T_T[t,])
                                          - J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) - a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G)
                                          - t(J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G))
                                          + new_G %*% (P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])) %*% t(new_G)))
      
    }else{
      new_W = parameters_current$W
    }
    
    parameters_current$G = new_G
    parameters_current$W = new_W
    
    if(log[[1]])
    {
      log_likelihood_history_subset = log_likelihood_history[1:iteration]
      G_history_subset     = G_history[1:iteration, , ]
      W_history_subset     = W_history[1:iteration, , ]
      Sigma_history_subset = Sigma_history[1:iteration, , ]
      
      EM_subset_output = list( log_likelihood_history = log_likelihood_history_subset, 
                               G_history = G_history_subset, 
                               W_history = W_history_subset,
                               Sigma_history = Sigma_history_subset)
      
      save(EM_subset_output, file = paste(paste(log[[2]][1], sep = "", log[[2]][2]), sep = "", "_subset_EM.Rdata"))
      
    }
    
  }
  
  print(paste("END"))
  filter_output = filter(data, parameters_current)
  smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
  log_likelihood_history[iteration+1] = filter_output$log_likelihood
  print(filter_output$log_likelihood)
  
  if(learn_G)
  {
    print("Matrix G")
    print(parameters_current$G)
    G_history[iteration+1,,] = parameters_current$G
  }
  
  if(learn_W)
  {
    print("Matrix W")
    print(parameters_current$W)
    W_history[iteration+1,,] = parameters_current$W
  }
  
  if(learn_Sigma)
  {
    print("Matrix Sigma")
    print(parameters_current$Sigma)
    Sigma_history[iteration+1,,] = parameters_current$Sigma
  }
  
  return(list(log_likelihood_history = log_likelihood_history, 
              G_history = G_history, 
              W_history = W_history,
              Sigma_history = Sigma_history))
}

