rm(list=ls())
# Set the directory to the path of the current script

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path)) {
  setwd(dirname(script_path))
  message("Working directory set to: ", getwd())
} else {
  message("Script path not detected. Working directory not changed.")
}

source("../scripts/utils.R")
source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_filter_SMC.R")

source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../scripts/multinomial_probit_EM_Sigma.R")

library(tidyr)
library(dplyr)
library(gridExtra)
library(transport)
library(expm) # for sqrtm (matrix square root)
library(future.apply)

# SET UP THE SIMULATION
L = 3 # number of choices
Time= 10 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 15

name_init="trace_constrained_repeat_EM_p_"

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(678)

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

G = 0.5*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.25
  G[i-1,i] = sign*0.25
}

W = 0.05*diag(beta_size)

Sigma = 0.2*diag(L)
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.1
  Sigma[i-1,i] = sign*0.1
}

invSigma = solve(Sigma)

target_trace = sum(diag(invSigma))

# GENERATE DATA

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# EM algorithm
EM_iterations = 200

n_samples = 50

Sigma_repeated_sample  = array(NA, dim = c(n_samples, L, L))
log_likelihood_history = array(NA, dim = c(n_samples, EM_iterations+1))
log_likelihood_true    = array(NA, dim = c(n_samples))

# Filter
EP_max_iterations = 5
MC_samples = 10000
eps_sample_size = 10000
tol = 10e-3
EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}

for(i in 1:n_samples)
{
  print(i)
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = parameters$G, 
                            W = parameters$W, 
                            Sigma = diag(L))
  
  data_multinomial_probit_next = simulator_given_covariates(data_multinomial_probit, parameters, hyperparameters, Time)
  
  EM_output_EP = multinomial_probit_EM_invSigma(data_multinomial_probit_next, 
                                                EP_sampling_filter_to_use, type = "nonstd",
                                                initial_parameters, target_trace,
                                                EM_iterations, eps_sample_size, 
                                                log = c(FALSE, NA))
  
  EP_filter_output_truth   = EP_sampling_filter_to_use(data_multinomial_probit_next, parameters)
  
  Sigma_repeated_sample[i,,] = EM_output_EP$Sigma_history[EM_iterations+1,,]
  log_likelihood_history[i,] = EM_output_EP$log_likelihood_history
  log_likelihood_true[i] = EP_filter_output_truth$log_likelihood 
  
}

save(Sigma_repeated_sample,  file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_Sigma_repeated_sample.Rdata"))
save(log_likelihood_history, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_history.Rdata"))
save(log_likelihood_true,    file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_true.Rdata"))