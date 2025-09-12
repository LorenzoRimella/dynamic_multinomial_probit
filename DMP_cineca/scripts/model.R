kernel_parameters <- function(beta, parameters)
{
  return(list(a_t = parameters$G %*% beta, P_t = parameters$W))
}

simulator_given_covariates <- function(data, parameters, hyperparameters = list(lambda_M, N_max, prediction_bool = FALSE), Time)
{
  Beta = matrix(NA, ncol = dim(parameters$a_0)[1], nrow =  Time+1)
  X    = list()
  Y    = list()
  prediction_prob = list()
  
  Beta[1,] = mvrnorm(n = 1, parameters$a_0, parameters$P_0)
  for(t in 1:Time){
    
    kernel_param_list = kernel_parameters(Beta[t,], parameters)
    
    Beta[t+1,] = mvrnorm(n = 1, kernel_param_list$a_t, kernel_param_list$P_t)
    
    X[[t]] = data$X[[t]]
    
    Xbeta = lapply(seq(1, dim(X[[t]])[1]), function(i) t(X[[t]][i,,]) %*% t(t(Beta[t+1,])))
    fixed_utilities = t(do.call(cbind, Xbeta))
    
    epsilon = mvrnorm(n = dim(fixed_utilities)[1], rep(0, dim(fixed_utilities)[2]), parameters$Sigma)
    utilities = fixed_utilities + epsilon
    utilities_0 = cbind(rep(0, dim(utilities)[1]), utilities)
    
    Y[[t]] = apply(utilities_0, 1, which.max)-1
    
    if(hyperparameters$prediction_bool){
      pred_mean   = lapply(seq(1, dim(X[[t]])[1]), function(i) lapply(seq(1, L), function(j) X_minus_l(X[[t]][i,,], j, L) %*% t(t(Beta[t+1,]))))
      pred_cov  = lapply(seq(1, dim(fixed_utilities)[1]), 
                         function(i) lapply(seq(1, dim(fixed_utilities)[2]), 
                                            function(j) V_minus_l(j, dim(fixed_utilities)[2]) %*% parameters$Sigma %*% t(V_minus_l(j, dim(fixed_utilities)[2]))))
      
      prediction_prob_t = sapply(seq(1, dim(fixed_utilities)[1]), function(i) sapply(seq(1, dim(fixed_utilities)[2]), 
                                                                                     function(j)  Phi(pred_mean[[i]][[j]], 
                                                                                                      pred_cov[[i]][[j]])))
      prediction_prob[[t]] = rbind((1-apply(prediction_prob_t, 2, sum)), prediction_prob_t)
    }
  }
  
  if(hyperparameters$prediction_bool){
    return(list(Beta = Beta, Y = Y, X = data$X, prediction_prob = prediction_prob))
  }
  else{
    return(list(Beta = Beta, Y = Y, X = data$X))
  }
}

simulator <- function(X_cov, Y_cov, parameters, hyperparameters = list(lambda_M, N_max, prediction_bool = FALSE), Time)
{
  L = dim(parameters$Sigma)[1]
  M_t = rpois(Time, hyperparameters$lambda_M) # number of customers per time steps
  M_t[M_t==0]=2
  M_t[M_t==1]=2
  M_t = apply(rbind(M_t, rep(hyperparameters$N_max, T)), 2, min) #  ensure it is lower than the total nr of customers
  
  customers_index_set = lapply( M_t , function(x){ sample(seq(1, hyperparameters$N_max, by = 1), x) }) # customer per time step
  property_index_set  = lapply( M_t , function(m_t) sapply(seq(1, m_t), function(x) sample(seq(1, hyperparameters$Y_max, by = 1), L)))
  
  Beta = matrix(NA, ncol = dim(parameters$a_0)[1], nrow =  Time+1)
  X    = list()
  Y    = list()
  prediction_prob = list()
  
  Beta[1,] = mvrnorm(n = 1, parameters$a_0, parameters$P_0)
  for(t in 1:Time){
    
    kernel_param_list = kernel_parameters(Beta[t,], parameters)
    
    Beta[t+1,] = mvrnorm(n = 1, kernel_param_list$a_t, kernel_param_list$P_t)
    
    custumers_cov  = X_cov[customers_index_set[[t]],]
    properties_cov = lapply(seq(1, M_t[[t]]), function(i) Y_cov[property_index_set[[t]][,i],])

    X_list = lapply(seq(1, M_t[[t]]), function(i) t(cbind(t(abind(lapply(seq(1, L), 
                                                                   function(l) t(one_hot(l,L))[1,] %x% custumers_cov[i,]), along = 2)), 
                                                           properties_cov[[i]], 
                                                           diag(L))))
    X_list = lapply(seq(1, M_t[[t]]), function(i) array(X_list[[i]], dim = c(1,dim(X_list[[i]]))))
    
    X[[t]] = abind(X_list, along = 1)
    
    Xbeta = lapply(seq(1, M_t[[t]]), function(i) t(X[[t]][i,,]) %*% t(t(Beta[t+1,])))
    fixed_utilities = t(do.call(cbind, Xbeta))
    
    epsilon = mvrnorm(n = dim(fixed_utilities)[1], rep(0, dim(fixed_utilities)[2]), parameters$Sigma)
    utilities = fixed_utilities + epsilon
    utilities_0 = cbind(rep(0, dim(utilities)[1]), utilities)
    
    Y[[t]] = apply(utilities_0, 1, which.max)-1
    
    if(hyperparameters$prediction_bool){
      pred_mean   = lapply(seq(1, dim(X[[t]])[1]), function(i) lapply(seq(1, L), function(j) X_minus_l(X[[t]][i,,], j, L) %*% t(t(Beta[t+1,]))))
      pred_cov  = lapply(seq(1, dim(fixed_utilities)[1]), 
                         function(i) lapply(seq(1, dim(fixed_utilities)[2]), 
                                            function(j) V_minus_l(j, dim(fixed_utilities)[2]) %*% parameters$Sigma %*% t(V_minus_l(j, dim(fixed_utilities)[2]))))
      
      prediction_prob_t = sapply(seq(1, dim(fixed_utilities)[1]), function(i) sapply(seq(1, dim(fixed_utilities)[2]), 
                                                                                     function(j)  Phi(pred_mean[[i]][[j]], 
                                                                                                      pred_cov[[i]][[j]])))
      prediction_prob[[t]] = rbind((1-apply(prediction_prob_t, 2, sum)), prediction_prob_t)
    }
  }
  
  if(hyperparameters$prediction_bool){
    return(list(Beta = Beta, Y = Y, X = X, prediction_prob = prediction_prob))
  }
  else{
    return(list(Beta = Beta, Y = Y, X = X))
  }
}
