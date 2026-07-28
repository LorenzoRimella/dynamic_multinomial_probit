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
Time= 10 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 10

name_init="full_parallel_filt_smoot_lik_p_"

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

alpha = 0.5
G = alpha*diag(beta_size)

W = (1-alpha*alpha)*P_0

Sigma = matrix(c(1.6, -0.1, 0,
                 -0.1, 0.2, 0.1,
                 0, 0.1, 1.2),
               nrow = L, byrow = TRUE)

# GENERATE DATA

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# COMPUTE FILTERS
EP_max_iterations = 11
tol= 1e-3
MC_samples = 1e4

sliding_window = Time + 1

# FIXED COVARIATES; RESIMULATED STATES AND RESPONSES ACROSS REPLICATES

library(future)
library(future.apply)

# Detect CPUs from SLURM environment
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Fall back to available cores if SLURM variable not set
if (is.na(slurm_cpus) || slurm_cpus < 1) {
  slurm_cpus <- parallel::detectCores(logical = FALSE)
}

message("Using ", slurm_cpus, " workers for parallel processing.")

# Use multisession (safer on HPC clusters)
plan(multisession, workers = slurm_cpus)

n_replicates = 100

results <- future_lapply(1:n_replicates, function(repl) {
  set.seed(12345+repl)
  print(repl)

  data_multinomial_probit_repl = simulator_given_covariates(data_multinomial_probit, parameters, hyperparameters, Time)

  EP_filter_sampler_output   = EP_sampling_filter(data_multinomial_probit_repl, parameters, EP_max_iterations, MC_samples, tol)
  EP_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit_repl, parameters, EP_filter_sampler_output)

  ADF_filter_sampler_output   = ADF_sampling_filter(data_multinomial_probit_repl, parameters, MC_samples)
  ADF_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit_repl, parameters, ADF_filter_sampler_output)

  exact_filter_sampler_output   = exact_sampling_filter(data_multinomial_probit_repl, parameters, sliding_window, MC_samples)
  exact_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit_repl, parameters, exact_filter_sampler_output)

  SMC_filter_sampler_output   = SMC_sampling_filter(data_multinomial_probit_repl, parameters, MC_samples)
  SMC_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit_repl, parameters, SMC_filter_sampler_output)

  df_EP = data.frame(
    replicate = repl,
    type = "EP",
    filter_running_time = EP_filter_sampler_output$running_time,
    smoother_running_time = EP_smoother_sampler_output$running_time,
    log_likelihood = EP_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = mean_wasserstein_sampling(EP_filter_sampler_output$sampling_history, exact_filter_sampler_output$sampling_history),
    smoothing_distance_from_exact = mean_wasserstein_gaussian(EP_smoother_sampler_output, exact_smoother_sampler_output)
  )

  df_ADF = data.frame(
    replicate = repl,
    type = "ADF",
    filter_running_time = ADF_filter_sampler_output$running_time,
    smoother_running_time = ADF_smoother_sampler_output$running_time,
    log_likelihood = ADF_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = mean_wasserstein_sampling(ADF_filter_sampler_output$sampling_history, exact_filter_sampler_output$sampling_history),
    smoothing_distance_from_exact = mean_wasserstein_gaussian(ADF_smoother_sampler_output, exact_smoother_sampler_output)
  )

  df_SMC = data.frame(
    replicate = repl,
    type = "SMC",
    filter_running_time = SMC_filter_sampler_output$running_time,
    smoother_running_time = SMC_smoother_sampler_output$running_time,
    log_likelihood = SMC_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = mean_wasserstein_sampling(SMC_filter_sampler_output$sampling_history, exact_filter_sampler_output$sampling_history),
    smoothing_distance_from_exact = mean_wasserstein_gaussian(SMC_smoother_sampler_output, exact_smoother_sampler_output)
  )

  df_exact = data.frame(
    replicate = repl,
    type = "exact",
    filter_running_time = exact_filter_sampler_output$running_time,
    smoother_running_time = exact_smoother_sampler_output$running_time,
    log_likelihood = exact_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = 0,
    smoothing_distance_from_exact = 0
  )

  return(list(df_EP, df_ADF, df_SMC, df_exact))
}, future.seed=TRUE)

data_frame_output <- do.call(rbind, unlist(results, recursive = FALSE))
save(data_frame_output, file = paste0(directory_name, name_simulation, "_data.Rdata"))

save(data_multinomial_probit, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_data_multinomial_probit_single.Rdata"))
