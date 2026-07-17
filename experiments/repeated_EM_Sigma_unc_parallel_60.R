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
library(transport)
library(expm) # for sqrtm (matrix square root)

library(future)
library(future.apply)

library(dplyr)
library(tidyr)


# SET UP THE SIMULATION
L = 3 # number of choices
Time= 60 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 15

name_init="parallel_Sigma_trace_AGL_repeat_EM_p_"

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

alpha = 0.5
G = alpha*diag(beta_size)

W = (1-alpha*alpha)*P_0

Sigma = diag(L)
Sigma[1,1] = 1.6
Sigma[2,2] = 0.2
Sigma[3,3] = 1.2
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.1
  Sigma[i-1,i] = sign*0.1
}

# GENERATE DATA

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# EM algorithm
EM_iterations = 200

n_samples = 48

Sigma_repeated_sample  = array(NA, dim = c(n_samples, L, L))
log_likelihood_history = array(NA, dim = c(n_samples, EM_iterations+1))
log_likelihood_true    = array(NA, dim = c(n_samples))

# Filter
EP_max_iterations = 5
MC_samples = 10000
eps_sample_size = 10000
tol = 10e-3
target_trace = L

EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}

# Detect CPUs from SLURM environment
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Fall back to available cores if SLURM variable not set
if (is.na(slurm_cpus) || slurm_cpus < 1) {
  slurm_cpus <- parallel::detectCores(logical = FALSE)
}

message("Using ", slurm_cpus, " workers for parallel processing.")

# Use multisession (safer on HPC clusters)
plan(multisession, workers = slurm_cpus)

results <- future_lapply(1:n_samples, function(i) {
  
  set.seed(678+i)
  print(i)
  
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = parameters$G, 
                            W = parameters$W, 
                            Sigma = diag(L))
  
  data_multinomial_probit_next = simulator_given_covariates(data_multinomial_probit, parameters, hyperparameters, Time)
  
  EM_output_EP = multinomial_probit_EM_Sigma(data_multinomial_probit_next, 
                                             EP_sampling_filter_to_use, type = "nonstd",
                                             initial_parameters, target_trace,
                                             EM_iterations, eps_sample_size, 
                                             log = c(FALSE, NA))
  
  EP_filter_output_truth   = EP_sampling_filter_to_use(data_multinomial_probit_next, parameters)
  
  return(list(EM_output_EP$Sigma_history[EM_iterations+1,,], 
              EM_output_EP$log_likelihood_history, 
              EP_filter_output_truth$log_likelihood))
  
}, future.seed=TRUE)

for(i in 1:n_samples)
{
  Sigma_repeated_sample[i,,] = results[[i]][[1]]
  log_likelihood_history[i,] = results[[i]][[2]]
  log_likelihood_true[i] = results[[i]][[3]]
}

save(Sigma_repeated_sample,  file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_Sigma_repeated_sample.Rdata"))
save(log_likelihood_history, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_history.Rdata"))
save(log_likelihood_true,    file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_true.Rdata"))
