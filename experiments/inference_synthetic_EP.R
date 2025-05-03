rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../scripts/utils.R")
source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_filter_SMC.R")
source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
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
Time= 5 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 1 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 20

name_init="inference_synthetic_EP_p_"

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(15)

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.5*(1/beta_size)*diag(beta_size)

G = 0.3*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.2
  G[i-1,i] = sign*0.2
}

W = 0.01*diag(beta_size)
is_positive_definite(W)

Sigma = 1.25*diag(L)
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.5
  Sigma[i-1,i] = sign*0.5
}
is_positive_definite(Sigma)

# GENERATE DATA
# Trace=1
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

# PLOT DATA
# Convert Beta matrix to long format
df_beta      = as.data.frame(data_multinomial_probit$Beta)
df_beta$Time = seq_len(Time+1)  

df_beta_long = df_beta %>%
  gather(Column, Value, -Time) %>%
  mutate(Type = "Beta")  

time_steps <- seq_along(data_multinomial_probit$Y)  
data_list <- lapply(seq_along(data_multinomial_probit$Y), 
                    function(i) {data.frame(table(data_multinomial_probit$Y[[i]])) %>%
                        mutate(time = i)})

df_y <- do.call(rbind, data_list)
names(df_y) <- c("Category", "Count", "Time")  
df_y$Time <- as.factor(df_y$Time)
df_y$Category <- factor(df_y$Category, levels = as.character(sort(as.numeric(levels(df_y$Category)))))

plot_1 = ggplot(df_y, aes(x = Time, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  labs(x = "Time", y = "Count", fill = "Category") +
  theme_minimal() 

plot_2 =ggplot(df_beta_long, aes(x = Time, y = Value, color = Type)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "Time", y = "") +
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  theme_minimal() 

grid.arrange(plot_1, plot_2, ncol = 1)

# True likelihood
EP_max_iterations = 1000
MC_samples = 10000
tol = 10e-3

EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}

EP_filter_output_truth   = EP_sampling_filter_to_use(data_multinomial_probit, parameters)
EP_smoother_output_truth = Gaussian_smoother(data_multinomial_probit, parameters, EP_filter_output_truth)

# EM algorithm
EM_iterations = 2000
learn_G = FALSE
learn_W = FALSE
learn_Sigma = TRUE

if((learn_G) & !(learn_W) & !(learn_Sigma))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = parameters$W, 
                            Sigma = parameters$Sigma)
}


if((learn_G & learn_W) & !(learn_Sigma))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = diag(dim(parameters$W)[1]), 
                            Sigma = parameters$Sigma)
}

if(learn_Sigma & !(learn_G) & !(learn_W))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = parameters$G, 
                            W = parameters$W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & !(learn_W))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = parameters$W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & learn_W)
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = diag(dim(parameters$W)[1]), 
                            Sigma = diag(L))
}

EM_output_EP = multinomial_probit_EM(data_multinomial_probit, EP_sampling_filter_to_use, initial_parameters,
                                  EM_iterations,
                                  learn_G, learn_W, learn_Sigma, 
                                  type = "nonstd",
                                  log = list(TRUE, c(directory_name, name_simulation)))

save(data_multinomial_probit, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_nonstd_balanced_data_multinomial_probit.Rdata"))
save(EP_filter_output_truth, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_nonstd_balanced_EP_truth.Rdata"))
save(EM_output_EP,           file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_nonstd_balanced_EP_EM.Rdata"))

Trace_EM = sum(diag(EM_output_EP$Sigma_history[EM_iterations+1,,]))
EM_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0/Trace_EM, 
                     G = parameters$G, 
                     W = parameters$W/Trace_EM, 
                     Sigma = EM_output_EP$Sigma_history[EM_iterations+1,,]/Trace_EM)

EP_filter_output_EM   = EP_sampling_filter_to_use(data_multinomial_probit, EM_parameters)
EP_smoother_output_EM = Gaussian_smoother(data_multinomial_probit, EM_parameters, EP_filter_output_EM)

# PLOT LOG-LIKELIHOOD
# Create a data frame for the line plot
df <- data.frame(
  iterations = 1:(EM_iterations+1),
  log_likelihood = EM_output_EP$log_likelihood_history
)

# Create the plot
ggplot() +
  geom_line(data = df, aes(x = iterations, y = log_likelihood, color = "EM Log-Likelihood"), size = 1) + # Line plot
  geom_hline(aes(yintercept = EP_filter_output_truth$log_likelihood, color = "EP Log-Likelihood"), 
             linetype = "dashed", linewidth = 1) + # Horizontal line
  labs(
    x = "Iterations",  # X-axis label
    y = "Log-Likelihood",  # Y-axis label
    color = "Legend"  # Legend title
  ) +
  scale_color_manual(
    values = c("EM Log-Likelihood" = "blue", "EP Log-Likelihood" = "red"), # Custom colors
    labels = c("EM Log-Likelihood" = "EM Log-Likelihood", "EP Log-Likelihood" = "EP Exact Output") # Custom labels
  ) +
  theme_minimal()  # Minimal theme
  
########################################################################
####################################NOT BALANCED
########################################################################
########################################################################

set.seed(5)

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.5*(1/beta_size)*diag(beta_size)

G = 0.3*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.2
  G[i-1,i] = sign*0.2
}

W = 0.01*diag(beta_size)
is_positive_definite(W)

Sigma = 1.25*diag(L)
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.5
  Sigma[i-1,i] = sign*0.5
}
is_positive_definite(Sigma)

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


# GENERATE DATA
# Trace=1
parameters = list(a_0 = a_0, P_0 = L*P_0/sum(diag(Sigma)), G = G, W = L*W/sum(diag(Sigma)), Sigma = L*Sigma/sum(diag(Sigma)))

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# PLOT DATA
# Convert Beta matrix to long format
df_beta      = as.data.frame(data_multinomial_probit$Beta)
df_beta$Time = seq_len(Time+1)  

df_beta_long = df_beta %>%
  gather(Column, Value, -Time) %>%
  mutate(Type = "Beta")  

time_steps <- seq_along(data_multinomial_probit$Y)  
data_list <- lapply(seq_along(data_multinomial_probit$Y), 
                    function(i) {data.frame(table(data_multinomial_probit$Y[[i]])) %>%
                        mutate(time = i)})

df_y <- do.call(rbind, data_list)
names(df_y) <- c("Category", "Count", "Time")  
df_y$Time <- as.factor(df_y$Time)
df_y$Category <- factor(df_y$Category, levels = as.character(sort(as.numeric(levels(df_y$Category)))))

plot_1 = ggplot(df_y, aes(x = Time, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  labs(x = "Time", y = "Count", fill = "Category") +
  theme_minimal() 

plot_2 =ggplot(df_beta_long, aes(x = Time, y = Value, color = Type)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "Time", y = "") +
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  theme_minimal() 

grid.arrange(plot_1, plot_2, ncol = 1)

# True likelihood

EP_max_iterations = 1000
MC_samples = 10000
tol = 10e-3

EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}


EP_filter_output_truth   = EP_sampling_filter_to_use(data_multinomial_probit, parameters)
EP_smoother_output_truth = Gaussian_smoother(data_multinomial_probit, parameters, EP_filter_output_truth)

# EM algorithm
EM_iterations = 2000
learn_G = FALSE
learn_W = FALSE
learn_Sigma = TRUE

if((learn_G) & !(learn_W) & !(learn_Sigma))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = parameters$W, 
                            Sigma = parameters$Sigma)
}


if((learn_G & learn_W) & !(learn_Sigma))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = diag(dim(parameters$W)[1]), 
                            Sigma = parameters$Sigma)
}

if(learn_Sigma & !(learn_G) & !(learn_W))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = parameters$G, 
                            W = parameters$W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & !(learn_W))
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = parameters$W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & learn_W)
{
  initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = diag(beta_size), 
                            W = diag(dim(parameters$W)[1]), 
                            Sigma = diag(L))
}

EM_output_EP = multinomial_probit_EM(data_multinomial_probit, EP_sampling_filter_to_use, initial_parameters,
                                     EM_iterations,
                                     learn_G, learn_W, learn_Sigma, 
                                     type = "nonstd",
                                     log = list(TRUE, c(directory_name, name_simulation)))

save(data_multinomial_probit, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_nonstd_unbalanced_2_data_multinomial_probit.Rdata"))
save(EP_filter_output_truth, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_nonstd_unbalanced_2_EP_truth.Rdata"))
save(EM_output_EP,           file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_nonstd_unbalanced_2_EP_EM.Rdata"))


Trace_EM = sum(diag(EM_output_EP$Sigma_history[EM_iterations+1,,]))
EM_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0/Trace_EM, 
                     G = parameters$G, 
                     W = parameters$W/Trace_EM, 
                     Sigma = EM_output_EP$Sigma_history[EM_iterations+1,,]/Trace_EM)

EP_filter_output_EM   = EP_sampling_filter_to_use(data_multinomial_probit, EM_parameters)
EP_smoother_output_EM = Gaussian_smoother(data_multinomial_probit, EM_parameters, EP_filter_output_EM)

# PLOT LOG-LIKELIHOOD
# Create a data frame for the line plot
df <- data.frame(
  iterations = 1:(EM_iterations+1),
  log_likelihood = EM_output_EP$log_likelihood_history
)

# Create the plot
ggplot() +
  geom_line(data = df, aes(x = iterations, y = log_likelihood, color = "EM Log-Likelihood"), size = 1) + # Line plot
  geom_hline(aes(yintercept = EP_filter_output_truth$log_likelihood, color = "EP Log-Likelihood"), 
             linetype = "dashed", linewidth = 1) + # Horizontal line
  labs(
    x = "Iterations",  # X-axis label
    y = "Log-Likelihood",  # Y-axis label
    color = "Legend"  # Legend title
  ) +
  scale_color_manual(
    values = c("EM Log-Likelihood" = "blue", "EP Log-Likelihood" = "red"), # Custom colors
    labels = c("EM Log-Likelihood" = "EM Log-Likelihood", "EP Log-Likelihood" = "EP Exact Output") # Custom labels
  ) +
  theme_minimal()  # Minimal theme



###############################################################################
###############################################################################
###############################################################################
# Set up a 3x3 plot grid

load("~/R_directory/dynamic_multinomial_probit/data/inference_synthetic_EP_p_7_mt_20_L_3_Time_5/inference_synthetic_EP_p_7_mt_20_L_3_Time_5_nonstd_unbalanced_2_EP_EM.Rdata")
EM_output_EP_unbalanced = EM_output_EP

load("~/R_directory/dynamic_multinomial_probit/data/inference_synthetic_EP_p_7_mt_20_L_3_Time_5/inference_synthetic_EP_p_7_mt_20_L_3_Time_5_nonstd_balanced_EP_EM.Rdata")
EM_output_EP_balanced = EM_output_EP

par(mfrow = c(3, 3), mar = c(2, 2, 2, 1), oma = c(0, 0, 2, 0))  # optional: adjust margins

# Plot each matrix element
for (i in 1:3) {
  for (j in 1:3) {
    
    balanced   = EM_output_EP_balanced$Sigma_history[800:1000,i,j]
    unbalanced = EM_output_EP_unbalanced$Sigma_history[800:1000,i,j]
    
    unbalanced = abs(unbalanced - unbalanced[1])
    balanced   = abs(balanced   - balanced[1])
    
    # Define y-limits to fit both
    ylims <- range(c(balanced, unbalanced))
    
    # Plot
    plot(unbalanced, type = "l", col = "red", lwd = 2,
         ylim = ylims, main = paste0("(", i, ",", j, ")"),
         xlab = "", ylab = "")
    lines(balanced, col = "blue", lwd = 2)

    # Optional: add legend only in the last subplot
    if (i == 1 && j == 1) {
      legend("bottomright", legend = c("Unbalanced", "Balanced"),
             col = c("red", "blue"), lty = 1, lwd = 2, bty = "n", cex = 0.8)
    }
  }
}
    
    
    
    
    
    