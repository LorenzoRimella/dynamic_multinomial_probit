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

# SET UP THE SIMULATION
L = 3 # number of choices
Time= 50 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 50

name_init="filt_smoot_lik_traj_p_"

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(12345)

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

# GENERATE DATA

parameters = list(a_0 = a_0, P_0 = L*P_0/sum(diag(Sigma)), G = G, W = L*W/sum(diag(Sigma)), Sigma = L*Sigma/sum(diag(Sigma)))

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# COMPUTE FILTERS
EP_max_iterations = 11
tol= 1e-3
MC_samples = 1e4

sliding_window = 10

set.seed(12345)
  
EP_filter_sampler_output   = EP_sampling_filter(data_multinomial_probit, parameters, EP_max_iterations, MC_samples, tol)
EP_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, EP_filter_sampler_output)
  
ADF_filter_sampler_output   = ADF_sampling_filter(data_multinomial_probit, parameters, MC_samples)
ADF_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, ADF_filter_sampler_output)
  
save(data_multinomial_probit, file = paste(paste(directory_name, sep = "", 
                                                 name_simulation), sep = "", "_data_multinomial_probit_single.Rdata"))

save(EP_filter_sampler_output, file = paste(paste(directory_name, sep = "", 
                                                 name_simulation), sep = "", "_data_multinomial_probit_EP_filter.Rdata"))
save(EP_smoother_sampler_output, file = paste(paste(directory_name, sep = "", 
                                                 name_simulation), sep = "", "_data_multinomial_probit_EP_smoother.Rdata"))

save(ADF_filter_sampler_output, file = paste(paste(directory_name, sep = "", 
                                                 name_simulation), sep = "", "_data_multinomial_probit_ADF_filter.Rdata"))
save(ADF_smoother_sampler_output, file = paste(paste(directory_name, sep = "", 
                                                 name_simulation), sep = "", "_data_multinomial_probit_ADF_smoother.Rdata"))