library(dplyr)
library(tidyr)

library(ggplot2)


set.seed(678)

# SET UP THE SIMULATION
L = 3 # number of choices
Time= 10 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 15

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

alpha = 0.5
G = alpha*diag(beta_size)

W = (1-alpha*alpha)*P_0

Sigma = 0.2*diag(L)
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.1
  Sigma[i-1,i] = sign*0.1
}

# invSigma = solve(Sigma)

# target_trace = sum(diag(invSigma))

Sigma_time_mean = array(NA, c(3, 3, 3))
Sigma_time_sd = array(NA, c(3, 3, 3))
Sigma_time_quantile = array(NA, c(3, 2, 3, 3))

n_samples = 48

listing = c(10, 30, 60)

for(i in 1:3)
{
  Time = listing[i]
  string = paste("DMP_cineca/data/parallel_next1_Sigma_unconstrained_repeat_EM_p_9_mt_15_L_3_Time_", Time, sep = "")
  string = paste(string, "/parallel_next1_Sigma_unconstrained_repeat_EM_p_9_mt_15_L_3_Time_" ,sep = "")
  string = paste(string, Time ,sep = "")
  
  load(paste(string, "_log_likelihood_history.Rdata",sep = ""))
  load(paste(string, "_Sigma_repeated_sample.Rdata",sep = ""))
  load(paste(string, "_log_likelihood_true.Rdata",sep = ""))
  
  mean_Sigma = Reduce("+", lapply(seq(1, n_samples), function(i) Sigma_repeated_sample[i,,]))/n_samples
  sd_Sigma   = sqrt(Reduce("+", lapply(seq(1, n_samples), function(i) (Sigma_repeated_sample[i,,] - mean_Sigma)^2))/n_samples)
  
  Sigma_time_mean[i,,] =(mean_Sigma)
  Sigma_time_sd[i,,] = (sd_Sigma)
  
  Sigma_time_quantile[i,,,] = apply(Sigma_repeated_sample, c(2, 3), quantile, probs = c(0.025, 0.975))
  
  
}

EM_iterations = 200

# df <- as.data.frame(log_likelihood_history - t(t(log_likelihood_true))%*%array(1, dim = c(1, EM_iterations+1)))
# df$id <- factor(1:n_samples)  # Add an ID for each row

# Optional: give time-based column names
# colnames(df)[1:(EM_iterations+1)] <- paste0("t", 1:(EM_iterations+1))

# df_long <- pivot_longer(df, 
#                         cols = starts_with("t"), 
#                         names_to = "time", 
#                         values_to = "value")

# Convert 'time' to numeric (optional, for ordering)
# df_long$time <- as.numeric(sub("t", "", df_long$time))

# ggplot(df_long, aes(x = time, y = value, group = id, color = id)) +
#   geom_line() +
#   labs(x = "Time", y = "Value", title = "N Trajectories Over T Time Points") +
#   theme_minimal()

# Get row and column indices for each element
positions <- expand.grid(row = 1:L, col = 1:L)

# Flatten the matrices
vec1_10 <- as.vector(Sigma_time_mean[1,,])
vec2_10 <- as.vector(Sigma_time_sd[1,,])

vec1_30 <- as.vector(Sigma_time_mean[2,,])
vec2_30 <- as.vector(Sigma_time_sd[2,,])

vec1_60 <- as.vector(Sigma_time_mean[3,,])
vec2_60 <- as.vector(Sigma_time_sd[3,,])

vec3 <- as.vector(Sigma)

# Create the formatted table
result <- data.frame(
  element = paste(positions$row, positions$col, sep = ","),
  EM_10_mean_sd = paste(round(vec1_10,2), "+/-", round(vec2_10,2)),
  EM_30_mean_sd = paste(round(vec1_30,2), "+/-", round(vec2_30,2)),
  EM_60_mean_sd = paste(round(vec1_60,2), "+/-", round(vec2_60,2)),
  Sigma_value = round(vec3,2)
)

print(result)
# xtable(result)

# optional: filter only lower diagonal part of Sigma
indeces <- c(1,2,5,3,6,9)

mean_10 = vec1_10[indeces]
mean_30 = vec1_30[indeces]
mean_60 = vec1_60[indeces]

vec1_10 = as.vector(Sigma_time_quantile[1,1,,])
vec2_10 = as.vector(Sigma_time_quantile[1,2,,])

vec1_30 = as.vector(Sigma_time_quantile[2,1,,])
vec2_30 = as.vector(Sigma_time_quantile[2,2,,])

vec1_60 = as.vector(Sigma_time_quantile[3,1,,])
vec2_60 = as.vector(Sigma_time_quantile[3,2,,])

vec1_10 <- vec1_10[indeces]
vec2_10 <- vec2_10[indeces]

vec1_30 <- vec1_30[indeces]
vec2_30 <- vec2_30[indeces]

vec1_60 <- vec1_60[indeces]
vec2_60 <- vec2_60[indeces]

vec3 <- vec3[indeces]

# plot estimated values and real values
lb_10 <- vec1_10
ub_10 <- vec2_10

lb_30 <- vec1_30
ub_30 <- vec2_30

lb_60 <- vec1_60
ub_60 <- vec2_60

df <- data.frame(
  mean_10 = mean_10,
  lb_10 = lb_10,
  ub_10 = ub_10,
  mean_30 = mean_30,
  lb_30 = lb_30,
  ub_30 = ub_30,
  mean_60 = mean_60,
  lb_60 = lb_60,
  ub_60 = ub_60,
  trueSigma = vec3,
  id = 1:length(lb_10)  # add an ID to distinguish each interval
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

# --- Step 1: Reshape to long format ---
df_long <- df %>%
  pivot_longer(
    cols = matches("^(mean|lb|ub)_"),
    names_to = c(".value", "level"),
    names_pattern = "(mean|lb|ub)_(\\d+)"
  )

# --- Step 2: Prepare points (true mean + true Sigma) ---
points_df <- df_long %>%
  select(id, level, trueSigma, mean) %>%
  pivot_longer(
    cols = c(trueSigma, mean),
    names_to = "group",
    values_to = "x"
  ) %>%
  mutate(
    y = id + as.numeric(level)/100 - 0.2  # same vertical offset by level
  )

# --- Step 3: Plot ---
ggplot(df_long) +
  # Segments for intervals
  geom_segment(
    aes(
      x = lb, xend = ub,
      y = id + as.numeric(level)/100 - 0.2,
      yend = id + as.numeric(level)/100 - 0.2,
      color = level
    ),
    size = 1.5
  ) +
  
  # Points for estimated mean and true value
  geom_point(
    data = points_df,
    aes(x = x, y = y, color = group),
    size = 3
  ) +
  
  # Unified color legend (interval levels + mean + true value)
  scale_color_manual(
    name = "",
    values = c(
      "10" = "gray80",
      "30" = "gray60",
      "60" = "gray40",
      "mean" = "blue",
      "trueSigma" = "red"
    ),
    breaks = c("10", "30", "60", "mean", "trueSigma"),
    labels = c("T = 10", "T = 30", "T = 60",
               "Estimated mean", "True value")
  ) +
  
  # Styling
  theme_minimal() +
  labs(
    x = "Value",
    y = "",
    title = expression("True vs estimated values of " * Sigma)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_y_continuous(breaks = df$id, labels = y_labels)

