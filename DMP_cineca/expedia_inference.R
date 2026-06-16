rm(list=ls())
# Set the directory to the path of the current script

source("DMP_cineca/scripts/utils.R")
source("DMP_cineca/scripts/model.R")
source("DMP_cineca/scripts/multinomial_probit_filter_EP.R")
source("DMP_cineca/scripts/multinomial_probit_filter_SMC.R")
source("DMP_cineca/scripts/multinomial_probit_filter_ADF.R")
source("DMP_cineca/scripts/multinomial_probit_filter_exact.R")
source("DMP_cineca/scripts/multinomial_probit_gaussian_smoothing.R")
source("DMP_cineca/scripts/multinomial_probit_EM_Sigma.R")

library(tidyr)
library(dplyr)
library(transport)
library(expm) # for sqrtm (matrix square root)

library(future)
library(future.apply)

library(dplyr)
library(tidyr)


# LOAD DATA
data_X = read.csv("DMP_cineca/data/data_X10_rand.csv")
data_X$prop_location_score1 = as.double(scale(data_X$prop_location_score1))
data_X$log_price_usd = as.double(scale(data_X$log_price_usd))
data_X$prop_starrating = as.double(scale(data_X$prop_starrating))

data_Y = read.csv("DMP_cineca/data/data_Y10_rand.csv")

covariates_names = c( #"prop_brand_bool",   
                      "prop_location_score1",
                      #"promotion_flag",
                      "log_price_usd",
                      "prop_starrating")

L = 10
t_0  = min(data_X$t)
Time = max(data_X$t)

end   = 50 #Time - t_0 + 1

X = list()
Y = list()

set.seed(678)

percent_search = 1

current_t = t_0
for(i in 1:end)
{
  current_X = subset(data_X,    t == current_t)
  current_Y = subset(data_Y, t == current_t)
  searches  = unique(current_X$srch_id)
  
  n_searches = floor(percent_search*length(searches))
  
  if(n_searches<2){n_searches = 2}
  
  searches = sample(searches, n_searches, replace = FALSE, prob = NULL)
  
  X_array = array(NA, dim = c(length(searches), length(covariates_names)+L, L))
  Y_array = array(NA, dim = c(length(searches)))
  index_s=1
  for(s in searches)
  {
    #print(s)
    current_X_search = subset(current_X, srch_id == s)
    current_X_search = current_X_search[order(current_X_search$position), ]
    current_X_search = current_X_search[, covariates_names]
    
    current_Y_search = subset(current_Y, srch_id == s)
    current_Y_search = current_Y_search$booking_pos_sum
    
    X_array[index_s,,] = rbind(t(as.matrix(current_X_search)), diag(rep(1, L)))
    Y_array[index_s]   = current_Y_search
    index_s = index_s +1
  }
  
  X[[i]] = X_array
  Y[[i]] = Y_array
  
  current_t = current_t+1
}

# NAME SIMULATION
Time= length(X)

p = dim(X[[1]])[2]-L # number of subject dependent covariate

name_init="inference_expedia_EP_p_"

expected_m_t = mean(sapply(seq(1,length(X)), function(i) return(dim(X[[i]])[1])))

# CREATE AND PLOT DATA

data_multinomial_probit = list(X = X, Y = Y)

# Parameters
beta_size = p + L

# EM algorithm
EM_iterations = 100

loglikelihoods = array(NA, c(5, EM_iterations+1))
Sigmas         = array(NA, c(5, EM_iterations+1, 10, 10))

Time_data = 50
counter = 1
for(alpha_num in c(0.5, 0.75, 0.9, 0.95, 0.99)){
  name_init = paste("inference_expedia_alpha_", alpha_num, sep="")
  name_init = paste(name_init, "_EP_p_", sep="")
  
  expected_m_t = mean(sapply(seq(1,length(X)), function(i) return(dim(X[[i]])[1])))
  
  output_name = sim_name(name_init, L, Time, p, expected_m_t)
  name_simulation = output_name$name
  directory_name  = paste("DMP_cineca/data/", output_name$name_for_dir, sep ="")
  load(paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EM.Rdata"))
  
  loglikelihoods[counter,] = EM_output_EP$log_likelihood_history
  Sigmas[counter,,,] = EM_output_EP$Sigma_history
  
  counter = counter +1
  
}

plot(loglikelihoods[,EM_iterations+1], type = "l")

alpha_num = 0.99

name_init = paste("inference_expedia_longEM_alpha_", alpha_num, sep="")
name_init = paste(name_init, "_EP_p_", sep="")

expected_m_t = mean(sapply(seq(1,length(X)), function(i) return(dim(X[[i]])[1])))

output_name = sim_name(name_init, L, 50, p, 14.64)
name_simulation = output_name$name
directory_name  = paste("DMP_cineca/data/", output_name$name_for_dir, sep ="")
load(paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EM.Rdata"))

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

G = alpha_num*diag(beta_size)

W = (1-alpha_num*alpha_num)*P_0

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

Sigma_array <- EM_output_EP$Sigma_history
loglik <- EM_output_EP$log_likelihood_history[1:101]

# Convert array to a tidy data frame
df <- as.data.frame.table(Sigma_array, responseName = "value")
colnames(df) <- c("iteration", "i", "j", "value")

df <- df %>%
  mutate(
    iteration = as.integer(iteration),
    i = as.integer(i),
    j = as.integer(j)
  )

library(ggplot2)

ggplot(df, aes(x = iteration, y = value)) +
  geom_line(linewidth = 0.3) +
  facet_grid(i ~ j) +
  theme_minimal(base_size = 8) +
  labs(
    title = expression(Sigma ~ " history"),  # ← Greek Sigma
    x = "Iteration",
    y = expression(Sigma)
  ) +
  scale_x_continuous(breaks = c(0, 500)) +   # ← only show 0 and 500
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 6)
  )

log_df <- data.frame(
  iteration = seq_along(EM_output_EP$log_likelihood_history),
  log_likelihood = EM_output_EP$log_likelihood_history
)

plot1 =ggplot(log_df, aes(x = iteration, y = log_likelihood)) +
  geom_line() +
  theme_minimal(base_size = 5) +
  labs(
    x = "Iteration",
    y = "log-likelihood"
  ) +
  scale_x_continuous(breaks = c(0, 500)) 


# Extract the final 10x10 slice
Sigma_final <- EM_output_EP$Sigma_history[501, , ]
Sigma_cor = cov2cor(Sigma_final)

# Create grid of positions
positions <- expand.grid(row = 1:L, col = 1:L)

# Flatten the matrix into a vector (row-major order)
positions$value <- as.vector(Sigma_cor)

# Plot heatmap
plot2=ggplot(positions, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "grey80") +
  scale_fill_viridis_c(name = "Value") +
  coord_fixed() +
  scale_x_continuous(breaks = 1:L) +
  scale_y_reverse(breaks = 1:L) +  # reverse to match image() orientation
  labs(
    x = "Row",
    y = "Col"
  ) +
  theme_minimal(base_size = 5) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 5),
    plot.title = element_text(hjust = 0.5)
  )

library(patchwork)

plot1 | plot2


L = 10
t_0  = min(data_X$t)
Time = max(data_X$t)

end   = Time - t_0 + 1

X = list()
Y = list()

set.seed(678)

percent_search = 1

current_t = t_0
for(i in 1:end)
{
  current_X = subset(data_X,    t == current_t)
  current_Y = subset(data_Y, t == current_t)
  searches  = unique(current_X$srch_id)
  
  n_searches = floor(percent_search*length(searches))
  
  if(n_searches<2){n_searches = 2}
  
  searches = sample(searches, n_searches, replace = FALSE, prob = NULL)
  
  X_array = array(NA, dim = c(length(searches), length(covariates_names)+L, L))
  Y_array = array(NA, dim = c(length(searches)))
  index_s=1
  for(s in searches)
  {
    #print(s)
    current_X_search = subset(current_X, srch_id == s)
    current_X_search = current_X_search[order(current_X_search$position), ]
    current_X_search = current_X_search[, covariates_names]
    
    current_Y_search = subset(current_Y, srch_id == s)
    current_Y_search = current_Y_search$booking_pos_sum
    
    X_array[index_s,,] = rbind(t(as.matrix(current_X_search)), diag(rep(1, L)))
    Y_array[index_s]   = current_Y_search
    index_s = index_s +1
  }
  
  X[[i]] = X_array
  Y[[i]] = Y_array
  
  current_t = current_t+1
}

data_multinomial_probit = list(X = X, Y = Y)

# Filter
parameters = list(a_0 = a_0,
                  P_0 = P_0,
                  G   = G,
                  W   = W,
                  Sigma = EM_output_EP$Sigma_history[501, , ])

EP_filter_sampler_output = EP_sampling_filter_to_use(data_multinomial_probit, parameters)

a_hist <- EP_filter_sampler_output$a_history      # T x M matrix (means)
P_hist <- EP_filter_sampler_output$P_history      # T x M x M array (covariances)

T_ <- nrow(a_hist)
M_ <- ncol(a_hist)

# 1) Extract variances (diagonal) for each time step -> T x M matrix
vars_mat <- t(sapply(1:T_, function(t) sapply(1:M_, function(i) P_hist[t, i,i])))


# Optional: provide names for dimensions if missing
if (is.null(colnames(a_hist))) colnames(a_hist) <- paste0("dim", seq_len(M_))
colnames(vars_mat) <- colnames(a_hist)

# If you have real time stamps, set them here; otherwise use 1:T_
time_vec <- seq_len(T_)  # replace with your vector, e.g. time_vec <- EP_filter_sampler_output$time

# Build tidy frame with mean, variance, sd, and 95% CI
df <- as.data.frame(a_hist) |>
  mutate(.time = time_vec) |>
  pivot_longer(-.time, names_to = "dimension", values_to = "mean") |>
  left_join(
    as.data.frame(vars_mat) |>
      mutate(.time = time_vec) |>
      pivot_longer(-.time, names_to = "dimension", values_to = "var"),
    by = c(".time", "dimension")
  ) |>
  mutate(
    sd    = sqrt(pmax(var, 0)),           # guard tiny negative numeric noise
    lower = mean - 1.96 * sd,
    upper = mean + 1.96 * sd
  )

# Assuming you already have df built (with .time, dimension, mean, lower, upper, etc.)
first_three <- colnames(a_hist)[1:3]
remaining_dims <- setdiff(colnames(a_hist), first_three)

# ---- 1️⃣ Combined plot for the first 3 dimensions ----
df_first3 <- df %>%
  filter(dimension %in% first_three) %>%
  mutate(
    dimension = recode(dimension,
                       "dim1" = "prop_location_score1",
                       "dim2" = "log_price_usd",
                       "dim3" = "prop_starrating"
    )
  )

df_first3$.time <- as.Date(df_first3$.time, origin = "2013-03-01")

p_color <- ggplot(df_first3, aes(.time, mean, color = dimension, fill = dimension)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  labs(
    x = "",
    y = expression(beta[t]),
    color = "Dimension",
    fill = "Dimension"
  ) +
  theme_minimal(base_size = 12)

p_color  # show the combined plot


# ---- 2️⃣ Faceted plots for remaining dimensions ----
pos_names <- paste0("pos", seq_along(remaining_dims))
rename_map <- setNames(pos_names, remaining_dims)

df_rest <- df %>%
  filter(dimension %in% remaining_dims) %>%
  mutate(
    dimension = recode(dimension, !!!rename_map),
    dimension = factor(dimension, levels = pos_names)
  )


df_rest$.time <- as.Date(df_rest$.time, origin = "2013-03-01")

p_facet <- ggplot(df_rest, aes(.time, mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  facet_wrap(~ dimension, scales = "free_y", ncol = 5) +
  labs(x = "Date", y = expression(beta[t])) +
  theme_minimal(base_size = 12)

p_facet  # show the faceted plot

ggsave("expedia123.png", p_color, width = 12, height = 3, dpi = 300)
ggsave("expediarest.png", p_facet, width = 12, height = 6, dpi = 300)
