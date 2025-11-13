rm(list=ls())
# Set the directory to the path of the current script

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

covariates_names = c( "prop_brand_bool",   
                      "prop_location_score1",
                      "promotion_flag",
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

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

alpha = 0.5
G = alpha*diag(beta_size)

W = (1-alpha*alpha)*P_0

# the invariant distr is N(0,P_0)

Sigma = diag(L)

initial_parameters = list(a_0 = a_0,
                          P_0 = P_0,
                          G   = G,
                          W   = W,
                          Sigma = Sigma)

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

# EM algorithm
EM_iterations = 500

load("DMP_cineca/data/inference_expedia_EP_p_5_mt_14.64_L_10_Time_50/inference_expedia_EP_p_5_mt_14.64_L_10_Time_50_EM.Rdata")

Sigma_array <- EM_output_EP$Sigma_history
loglik <- EM_output_EP$log_likelihood_history

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
    y = "Value"
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

ggplot(log_df, aes(x = iteration, y = log_likelihood)) +
  geom_line() +
  theme_minimal(base_size = 10) +
  labs(
    title = expression("Log-likelihood history"),
    x = "Iteration",
    y = expression(log~L)
  ) +
  scale_x_continuous(breaks = c(0, 500)) 


# Extract the final 10x10 slice
Sigma_final <- EM_output_EP$Sigma_history[501, , ]

# Create grid of positions
positions <- expand.grid(row = 1:L, col = 1:L)

# Flatten the matrix into a vector (row-major order)
positions$value <- as.vector(Sigma_final)

# Plot heatmap
ggplot(positions, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "grey80") +
  scale_fill_viridis_c(name = "Value") +
  coord_fixed() +
  scale_x_continuous(breaks = 1:L) +
  scale_y_reverse(breaks = 1:L) +  # reverse to match image() orientation
  labs(
    title = expression(Sigma ~ " matrix"),
    x = "Row",
    y = "Col"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
