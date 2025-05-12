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

expected_m_t = 30

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
P_0 = diag(beta_size)

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

# parameters = list(a_0 = a_0, P_0 = L*P_0/sum(diag(Sigma)), G = G, W = L*W/sum(diag(Sigma)), Sigma = L*Sigma/sum(diag(Sigma)))
Sigma_orig <- Sigma
invSigma <- solve(Sigma)
invSigma <- L*invSigma/sum(diag(invSigma))
Sigma <- solve(invSigma)
parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

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
EM_iterations = 200
learn_G = FALSE
learn_W = FALSE
learn_Sigma = TRUE

n_samples = 1

Sigma_repeated_sample  = array(NA, dim = c(n_samples, L, L))
log_likelihood_history = array(NA, dim = c(n_samples, EM_iterations+1))
log_likelihood_true    = array(NA, dim = c(n_samples))


for(i in 1:n_samples){
  print(i)
  
  Sigma_init <- diag(L)
  Sigma_init <- Sigma_init + 0.7 - 0.7*diag(L)
  invSigma_init <- solve(Sigma_init)
  invSigma_init <- L*invSigma_init/sum(diag(invSigma_init))
  Sigma_init <- solve(invSigma_init)
  
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0,
                            G = parameters$G,
                            W = parameters$W,
                            Sigma = Sigma_init)
                            # Sigma = diag(L))

  data_multinomial_probit_next = simulator_given_covariates(data_multinomial_probit, parameters, hyperparameters, Time)
  
  res_true <- EP_sampling_filter(data_multinomial_probit_next, parameters)
  print(res_true$log_likelihood)

  EM_output_EP = multinomial_probit_EM_constrained(data_multinomial_probit_next, EP_sampling_filter, initial_parameters,
                                       EM_iterations,
                                       learn_G, learn_W, learn_Sigma, type = "nonstd",
                                       log = c(FALSE, NA),
                                       target_trace = L)

  EP_filter_output_truth   = EP_sampling_filter(data_multinomial_probit_next, parameters)

  # trace_Sigma = sum(diag(EM_output_EP$Sigma_history[EM_iterations+1,,]))
  # Sigma_repeated_sample[i,,] = L*EM_output_EP$Sigma_history[EM_iterations+1,,]/trace_Sigma
  Sigma_repeated_sample[i,,] = EM_output_EP$Sigma_history[EM_iterations+1,,]
  log_likelihood_history[i,] = EM_output_EP$log_likelihood_history
  log_likelihood_true[i] = EP_filter_output_truth$log_likelihood

}

# save(Sigma_repeated_sample,  file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_Sigma_repeated_sample.Rdata"))
# save(log_likelihood_history, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_history.Rdata"))
# save(log_likelihood_true,    file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_log_likelihood_true.Rdata"))

# 
# load("../data/repeat_EM_p_9_mt_15_L_3_Time_10/repeat_EM_p_9_mt_15_L_3_Time_10_log_likelihood_history.Rdata")
# load("../data/repeat_EM_p_9_mt_15_L_3_Time_10/repeat_EM_p_9_mt_15_L_3_Time_10_log_likelihood_true.Rdata")
# load("../data/repeat_EM_p_9_mt_15_L_3_Time_10/repeat_EM_p_9_mt_15_L_3_Time_10_Sigma_repeated_sample.Rdata")



df <- as.data.frame(log_likelihood_history - t(t(log_likelihood_true))%*%array(1, dim = c(1, EM_iterations+1)))
df$id <- factor(1:n_samples)  # Add an ID for each row

# Optional: give time-based column names
colnames(df)[1:(EM_iterations+1)] <- paste0("t", 1:(EM_iterations+1))

df_long <- pivot_longer(df, 
                        cols = starts_with("t"), 
                        names_to = "time", 
                        values_to = "value")

# Convert 'time' to numeric (optional, for ordering)
df_long$time <- as.numeric(sub("t", "", df_long$time))

ggplot(df_long, aes(x = time, y = value, group = id, color = id)) +
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

print(result)

# optional: filter only lower diagonal part of Sigma
indeces <- c(1,2,5,3,6,9)
vec1 <- vec1[indeces]
vec2 <- vec2[indeces]
vec3 <- vec3[indeces]

# plot estimated values and real values
lb <- vec1 - vec2
ub <- vec1 + vec2

df <- data.frame(
  lb = lb,
  ub = ub,
  mSigma = vec1,
  trueSigma = vec3,
  id = 1:length(lb)  # add an ID to distinguish each interval
)

# Create the labels as expressions like expression(Sigma[x[1]])
y_labels <- c(expression(Sigma[11]),
              expression(Sigma[21]),
              expression(Sigma[31]),
              expression(Sigma[12]),
              expression(Sigma[22]),
              expression(Sigma[32]),
              expression(Sigma[13]),
              expression(Sigma[23]),
              expression(Sigma[33]))

# optional
y_labels <- y_labels[indeces]

y_pos <- 1:length(vec1)
points_df <- data.frame(
  x = c(vec1, vec3),
  y = rep(y_pos, 2),
  group = factor(rep(c("Estimated mean", "True value"), each = length(y_pos)))
)

ggplot(df) +
  geom_segment(aes(x = lb, xend = ub, y = id, yend = id), size = 1.5, color = "gray") +
  geom_point(data = points_df, aes(x = x, y = y, color = group), size = 4) +
  
  # Legend & aesthetics
  scale_color_manual(values = c("Estimated mean" = "blue", "True value" = "red")) +
  theme_minimal() +
  labs(
    x = "Value",
    y = "",
    title = expression("True vs estimated values of " * Sigma),
    color = "Group"
  )+
  theme(
    plot.title = element_text(hjust = 0.5)
  )+
  scale_y_continuous(breaks = df$id, labels = y_labels)



# plot the two matrices
library(reshape2)

# For mat1
df1 <- melt(parameters$Sigma)
colnames(df1) <- c("row", "col", "value")
df1$matrix <- "True"

df2 <- melt(mean_Sigma)
colnames(df2) <- c("row", "col", "value")
df2$matrix <- "Estimated mean"

combined_df <- rbind(df1, df2)
combined_df$matrix <- factor(combined_df$matrix, levels = c("True", "Estimated mean"))

ggplot(combined_df, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_reverse() +
  facet_wrap(~ matrix) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = expression(Sigma*" true vs estimated heatmaps"),
    x = "Column",
    y = "Row",
    fill = "Value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
