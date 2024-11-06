library(MASS)
library(abind)

one_hot <- function(position, dimension)
{
  vector = matrix(0, nrow = dimension, ncol = 1)
  
  if(position<=dimension){
    vector[position] = 1 
  }
  
  return(vector)
}

sferical_transformation <- function(J, kappa, dimensions)
{
  psi = list()  
  B = list()
  D = list()
  for(i in 1:length(kappa))
  {
    supporting_kappa = sqrt(J[[i]])*cos(kappa[[i]])
    supporting_kappa[length(kappa[[i]])] = sqrt(J[[i]])*1
    
    psi_elem = supporting_kappa*cumprod(c(1, sin(kappa[[i]])[1:(length(kappa[[i]])-1)]))
    psi[[i]] = psi_elem
    
    B[[i]] = matrix(psi[[i]][1:(length(psi[[i]])-dimensions[[i]][1])], nrow = dimensions[[i]][1], ncol = dimensions[[i]][2], byrow = TRUE)
    D[[i]] = psi[[i]][(length(psi[[i]])-dimensions[[i]][1] + 1):length(psi[[i]])]
  }

  return(list(psi = psi, B = B, D = D)) 
}

kernel_parameters <- function(beta, parameters)
{
  return(list(a_t = parameters$G %*% beta, P_t = parameters$W))
}

simulator_v2 <- function(X_1, X_2, parameters, hyperparameters = list(lambda_M = 20, N_max = 40), Time)
{
  M_t = rpois(Time, lambda_M) # number of customers per time steps
  M_t = apply(rbind(M_t, rep(N_max, T)), 2, min) #  ensure it is lower than the total nr of customers
  
  customers_index_set = lapply(        M_t , function(x){ sample(seq(1, N_max, by = 1), x) }) # customer per time step
  choice_index_set    = lapply(rep(L, Time), function(x){ sample(seq(1, C_max, by = 1), x) }) # choices per time step
  
  Beta = matrix(NA, ncol = dim(parameters$a_0)[1], nrow =  Time+1)
  X    = list()
  Y    = list()
  Utilities = list()
  
  Beta[1,] = mvrnorm(n = 1, parameters$a_0, parameters$P_0)
  for(t in 1:Time){
    
    kernel_param_list = kernel_parameters(Beta[t,], parameters)
    
    Beta[t+1,] = mvrnorm(n = 1, kernel_param_list$a_t, kernel_param_list$P_t)
    
    X_1_t = X_1[customers_index_set[[t]],] # customers covariate at time t
    X_2_t = X_2[choice_index_set[[t]],] # choices covariate at time t 
    X_2_t = cbind(rep(1, dim(X_2_t)[1]), X_2_t)
    
    X_2_t_choices = t(sapply(seq(1, dim(X_2_t)[1]), function(i) one_hot(i,dim(X_2_t)[1]) %x% X_2_t[i,]))
    
    # X_choices = lapply(seq(1, dim(X_2_t_choices)[1]), function(i) cbind(rep(1, dim(X_1_t)[1]), 
    #                                                                   cbind(X_1_t, t(t(rep(1, dim(X_1_t)[1]))) %*% X_2_t_choices[i,])))
    
    X_choices = lapply(seq(1, dim(X_2_t_choices)[1]), function(i) cbind(X_1_t, t(t(rep(1, dim(X_1_t)[1]))) %*% X_2_t_choices[i,]))
    
    X[[t]] = abind(X_choices, along = 3)
    
    Xbeta = lapply(seq(1, dim(X_2_t_choices)[1]), function(i) X_choices[[i]] %*% t(t(Beta[t+1,])))
    fixed_utilities = do.call(cbind, Xbeta)
    
    
    epsilon = mvrnorm(n = dim(fixed_utilities)[1], rep(0, dim(fixed_utilities)[2]), parameters$Sigma)
    utilities = fixed_utilities + epsilon
    utilities_0 = cbind(rep(0, dim(utilities)[1]), utilities)
    
    Y[[t]] = apply(utilities_0, 1, which.max)-1
    Utilities[[t]] = utilities_0
  }
  
  return(list(Beta = Beta, Y = Y, X = X, Utilities = Utilities))
}

simulator <- function(X_cov, parameters, hyperparameters = list(lambda_M = 20, N_max = 40), Time)
{
  L = dim(parameters$Sigma)[1]
  M_t = rpois(Time, lambda_M) # number of customers per time steps
  M_t[M_t==0]=2
  M_t[M_t==1]=2
  M_t = apply(rbind(M_t, rep(N_max, T)), 2, min) #  ensure it is lower than the total nr of customers
  
  customers_index_set = lapply(        M_t , function(x){ sample(seq(1, N_max, by = 1), x) }) # customer per time step
  
  Beta = matrix(NA, ncol = dim(parameters$a_0)[1], nrow =  Time+1)
  X    = list()
  Y    = list()
  prediction_prob = list()
  
  Beta[1,] = mvrnorm(n = 1, parameters$a_0, parameters$P_0)
  for(t in 1:Time){
    
    kernel_param_list = kernel_parameters(Beta[t,], parameters)
    
    Beta[t+1,] = mvrnorm(n = 1, kernel_param_list$a_t, kernel_param_list$P_t)
    
    X_1_t = X_cov[customers_index_set[[t]],] # customers covariate at time t
    X_2_t_choices = diag(rep(1, L)) # choices covariate at time t 
    
    X_choices = lapply(seq(1, dim(X_2_t_choices)[1]), function(i) cbind(X_1_t, t(t(rep(1, dim(X_1_t)[1]))) %*% X_2_t_choices[i,]))
    
    X[[t]] = abind(X_choices, along = 3)
    
    Xbeta = lapply(seq(1, dim(X_2_t_choices)[1]), function(i) X_choices[[i]] %*% t(t(Beta[t+1,])))
    fixed_utilities = do.call(cbind, Xbeta)
    
    epsilon = mvrnorm(n = dim(fixed_utilities)[1], rep(0, dim(fixed_utilities)[2]), parameters$Sigma)
    utilities = fixed_utilities + epsilon
    utilities_0 = cbind(rep(0, dim(utilities)[1]), utilities)
    
    Y[[t]] = apply(utilities_0, 1, which.max)-1
    
    pred_mean   = lapply(seq(1, dim(X[[t]])[1]), function(i) lapply(seq(1, L), function(j) X_minus_l(X[[t]][i,,], j, L) %*% t(t(Beta[t+1,]))))
    pred_cov  = lapply(seq(1, dim(fixed_utilities)[1]), 
                       function(i) lapply(seq(1, dim(fixed_utilities)[2]), 
                                          function(j) V_minus_l(j, dim(fixed_utilities)[2]) %*% parameters$Sigma %*% t(V_minus_l(j, dim(fixed_utilities)[2]))))
    
    prediction_prob_t = sapply(seq(1, dim(fixed_utilities)[1]), function(i) sapply(seq(1, dim(fixed_utilities)[2]), 
                                                                                 function(j)  Phi(pred_mean[[i]][[j]], 
                                                                                                  pred_cov[[i]][[j]])))
    prediction_prob[[t]] = rbind((1-apply(prediction_prob_t, 2, sum)), prediction_prob_t)
    # mean_prediction_prob_t = apply(prediction_prob_t, 1, mean)
    # prediction_prob[t,] = mean_prediction_prob_t
  }
  
  return(list(Beta = Beta, Y = Y, X = X, prediction_prob = prediction_prob))
}

comp_utilities <- function(X, Beta)
{
  Utilities = list()
  for(t in 1:Time){
    
    Beta_t = Beta[t+1,]
    
    X_t = X[[t]]
    
    Xbeta = lapply(seq(1, dim(X_t)[3]), function(i) X_t[,,i] %*% t(t(Beta_t)))
    fixed_utilities = do.call(cbind, Xbeta)

    Utilities[[t]] = fixed_utilities
  }
  
  return(list(Utilities = Utilities))
}
