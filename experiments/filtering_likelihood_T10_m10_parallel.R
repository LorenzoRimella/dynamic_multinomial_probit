rm(list = ls())

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

# Simulation setting
L = 3
Time = 10

p_X = 2
p_Y = 3
p = p_X*L + p_Y

expected_m_t = 10

name_init = "filtering_comparison_p_"

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(12345)

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

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

EP_max_iterations = 11
tol = 1e-3
MC_samples = 1e4

sliding_window = Time + 1

library(future)
library(future.apply)

# Detect CPUs from SLURM environment
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Fall back to available cores if SLURM variable not set
if (is.na(slurm_cpus) || slurm_cpus < 1) {
  slurm_cpus <- parallel::detectCores(logical = FALSE)
}

message("Using ", slurm_cpus, " workers for parallel processing.")

plan(multisession, workers = slurm_cpus)

n_replicates = 48

results <- future_lapply(1:n_replicates, function(repl) {
  set.seed(12345+repl)
  message("Replication ", repl)

  data_multinomial_probit_repl = simulator_given_covariates(data_multinomial_probit, parameters, hyperparameters, Time)

  EP_filter_sampler_output   = EP_sampling_filter(data_multinomial_probit_repl, parameters, EP_max_iterations, MC_samples, tol)
  
  ADF_filter_sampler_output   = ADF_sampling_filter(data_multinomial_probit_repl, parameters, MC_samples)

  exact_filter_sampler_output   = exact_sampling_filter(data_multinomial_probit_repl, parameters, sliding_window, MC_samples)

  SMC_filter_sampler_output   = SMC_sampling_filter(data_multinomial_probit_repl, parameters, MC_samples)

  df_EP = data.frame(
    replicate = repl,
    type = "EP",
    filter_running_time = EP_filter_sampler_output$running_time,
    log_likelihood = EP_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = mean_wasserstein_sampling(EP_filter_sampler_output$sampling_history, exact_filter_sampler_output$sampling_history),
    other = mean(array(unlist(EP_filter_sampler_output$sweep_history)))
  )

  df_ADF = data.frame(
    replicate = repl,
    type = "ADF",
    filter_running_time = ADF_filter_sampler_output$running_time,
    log_likelihood = ADF_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = mean_wasserstein_sampling(ADF_filter_sampler_output$sampling_history, exact_filter_sampler_output$sampling_history),
    other = 0
  )

  df_SMC = data.frame(
    replicate = repl,
    type = "SMC",
    filter_running_time = SMC_filter_sampler_output$running_time,
    log_likelihood = SMC_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = mean_wasserstein_sampling(SMC_filter_sampler_output$sampling_history, exact_filter_sampler_output$sampling_history),
    other = mean(array(unlist(SMC_filter_sampler_output$ESS_history)))
  )

  df_exact = data.frame(
    replicate = repl,
    type = "exact",
    filter_running_time = exact_filter_sampler_output$running_time,
    log_likelihood = exact_filter_sampler_output$log_likelihood,
    filter_distance_from_exact = 0,
    other = 0 
  )

  return(list(df_EP, df_ADF, df_SMC, df_exact))
}, future.seed=TRUE)

data_frame_output <- do.call(rbind, unlist(results, recursive = FALSE))
save(data_frame_output, file = paste0(directory_name, name_simulation, "_data.Rdata"))
