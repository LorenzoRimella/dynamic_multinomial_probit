rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# check the values of P_0 etc. to make sure that there are not deterministic cases
# repeat the EM algorithm and report the covergence traceplots 

source("../scripts/utils.R")
source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_filter_SMC.R")

source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../scripts/multinomial_probit_EM.R")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(transport)
library(expm) # for sqrtm (matrix square root)
library(future.apply)

library(dplyr)
library(tidyr)


# SET UP THE SIMULATION
L = 3 # number of choices
Time= 10 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 15

name_init="repeat_EM_p_"

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
P_0 = 1/p*diag(beta_size)

W = 0.05*diag(beta_size)

G = 0.5*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.25
  G[i-1,i] = sign*0.25
}

Sigma = 0.2*diag(L)
Sigma[1,1] = 0.5
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

# Check mean variance of the fixed effect
Sigma_cumul = parameters$P_0
sum_var     = array(0, dim= c(L,L))
elem = 0 
for(t in 1:Time)
{
  Sigma_cumul = parameters$G%*%Sigma_cumul%*%t(parameters$G) + parameters$W
  
  elem = elem + dim(data_multinomial_probit$X[[t]])[1]
  
  list_variances = lapply(seq(1, dim(data_multinomial_probit$X[[t]])[1]), function(i) t(data_multinomial_probit$X[[t]][i,,])%*%Sigma_cumul%*%data_multinomial_probit$X[[t]][i,,])
  
  sum_var = sum_var + Reduce("+", list_variances)
}

df <- expand.grid(x = 1:nrow(sum_var/elem), y = 1:ncol(sum_var/elem))
df$z <- as.vector(sum_var/elem)  # Assign values from the matrix

ggplot(df, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c() + 
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal()

# EM algorithm
EM_iterations = 100
learn_G = FALSE
learn_W = FALSE
learn_Sigma = TRUE

n_samples = 50

Sigma_repeated_sample  = array(NA, dim = c(n_samples, L, L))
log_likelihood_history = array(NA, dim = c(n_samples, EM_iterations+1))
log_likelihood_true    = array(NA, dim = c(n_samples))

for(i in 1:n_samples)
{
  print(i)
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = parameters$G, 
                            W = parameters$W, 
                            Sigma = diag(L))
  
  data_multinomial_probit_next = simulator_given_covariates(data_multinomial_probit, parameters, hyperparameters, Time)
  
  EM_output_EP = multinomial_probit_EM(data_multinomial_probit_next, EP_sampling_filter, initial_parameters,
                                       EM_iterations,
                                       learn_G, learn_W, learn_Sigma, type = "nonstd",
                                       log = c(FALSE, NA))
  
  EP_filter_output_truth   = EP_sampling_filter(data_multinomial_probit_next, parameters)
  
  # trace_Sigma = sum(diag(EM_output_EP$Sigma_history[EM_iterations+1,,]))
  # Sigma_repeated_sample[i,,] = L*EM_output_EP$Sigma_history[EM_iterations+1,,]/trace_Sigma
  Sigma_repeated_sample[i,,] = EM_output_EP$Sigma_history[EM_iterations+1,,]
  log_likelihood_history[i,] = EM_output_EP$log_likelihood_history
  log_likelihood_true[i] = EP_filter_output_truth$log_likelihood 
  
}

save(Sigma_repeated_sample,  file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_Sigma_repeated_sample.Rdata"))
save(log_likelihood_history, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_history.Rdata"))
save(log_likelihood_true,    file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_true.Rdata"))

df <- as.data.frame(log_likelihood_history - t(t(log_likelihood_true))%*%array(1, dim = c(1, EM_iterations+1)))
df$id <- factor(1:n_samples)  # Add an ID for each row

# Optional: give time-based column names
colnames(df)[1:(EM_iterations+1)] <- paste0("t", 1:(EM_iterations+1))

df_long <- pivot_longer(df, 
                        cols = starts_with("t"), 
                        names_to = "iteration", 
                        values_to = "diff")

# Convert 'time' to numeric (optional, for ordering)
df_long$iteration <- as.numeric(sub("t", "", df_long$iteration))

ggplot(df_long, aes(x = iteration, y = diff, group = id, color = 3)) +
  geom_line() +
  labs(x = "Time", y = "Value", title = "N Trajectories Over T Time Points") +
  theme_minimal()

matrix_normalization <- function(Sigma)
{
  trace_Sigma = sum(diag(Sigma))
  
  return(L*Sigma/trace_Sigma)
}

mean_Sigma = Reduce("+", lapply(seq(1, n_samples), function(i) matrix_normalization(Sigma_repeated_sample[i,,])))/n_samples
sd_Sigma   = sqrt(Reduce("+", lapply(seq(1, n_samples), function(i) (matrix_normalization(Sigma_repeated_sample[i,,]) - mean_Sigma)^2))/n_samples)

parameters$Sigma

# Get row and column indices for each element
positions <- expand.grid(row = 1:L, col = 1:L)

# Flatten the matrices
vec1 <- as.vector(mean_Sigma)
vec2 <- as.vector(sd_Sigma)
vec3 <- as.vector(parameters$Sigma)

# Create the formatted table
result <- data.frame(
  element = paste(positions$row, positions$col, sep = ","),
  EM_mean_sd = paste(round(vec1,2), "+/-", round(vec2,2)),
  Sigma_value = round(vec3,2)
)




