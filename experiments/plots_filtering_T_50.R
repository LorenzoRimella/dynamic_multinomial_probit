library(dplyr)
library(ggplot2)

load("data/filt_smoot_lik_p_9_mt_50_L_3_Time_50/filt_smoot_lik_p_9_mt_50_L_3_Time_50_data.Rdata")

# Compute mean and standard deviation for each type
summary_table <- data_frame_output %>%
  group_by(type) %>%
  summarise(across(
    where(is.numeric),
    list(mean = ~mean(.), sd = ~sd(.)),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(across(
    ends_with("_mean"),
    ~paste0(round(.x, 3), " ± ", round(get(sub("_mean$", "_sd", cur_column())), 3)),
    .names = "{.col}_formatted"
  )) %>%
  select(type, ends_with("_formatted")) %>%
  rename_with(~sub("_mean_formatted$", "", .x))

# Print table
print(summary_table, row.names = FALSE)

# plot trajectories
load("data/filt_smoot_lik_p_9_mt_50_L_3_Time_50/filt_smoot_lik_p_9_mt_50_L_3_Time_50_data_multinomial_probit.Rdata")
load("data/filt_smoot_lik_p_9_mt_50_L_3_Time_50/filt_smoot_lik_p_9_mt_50_L_3_Time_50_ADF_filter.Rdata")
load("data/filt_smoot_lik_p_9_mt_50_L_3_Time_50/filt_smoot_lik_p_9_mt_50_L_3_Time_50_EP_filter.Rdata")
load("data/filt_smoot_lik_p_9_mt_50_L_3_Time_50/filt_smoot_lik_p_9_mt_50_L_3_Time_50_SMC_filter.Rdata")
load("data/filt_smoot_lik_p_9_mt_50_L_3_Time_50/filt_smoot_lik_p_9_mt_50_L_3_Time_50_exact_filter.Rdata")

dim(data_multinomial_probit$X[[1]]) # n_obs x p x L
t(data_multinomial_probit$X[[1]][1,,]) # design matrix for first individual at first time

Time <- length(data_multinomial_probit$Y)

color_list <- c("blue", "black", "red", "magenta", "brown")

# individual-specific characteristics
indeces <- 1:6
# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
for (j in indeces) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = EP_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = EP_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j+3,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j+3,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = SMC_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = SMC_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(SMC_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = SMC_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(SMC_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+2],
    facet_row = j+3,
    facet_col = "SMC"
  )
  data_list[[length(data_list) + 1]] <- df_filter
}

plot_data <- bind_rows(data_list)

# Create the ggplot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), alpha = 0.5) +
  geom_line(aes(y = y2, color = color), alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(facet_row ~ facet_col) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    strip.text = element_text(size = 18),         # Facet labels (e.g., column/row labels)
    axis.title = element_text(size = 18),         # Axis titles (x and y)
    axis.text = element_text(size = 14)
  ) +
  labs(
    title = "Parameters for individual-specific characteristics",
    x = "Time",
    y = "Beta_t"
  )

# choice-specific characteristics
indeces <- 7:9
# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
for (j in indeces) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = EP_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = EP_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j+3,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j+3,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = SMC_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = SMC_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(SMC_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = SMC_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(SMC_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+2],
    facet_row = j+3,
    facet_col = "SMC"
  )
  data_list[[length(data_list) + 1]] <- df_filter
}

plot_data <- bind_rows(data_list)

# Create the ggplot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), alpha = 0.5) +
  geom_line(aes(y = y2, color = color), alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(facet_row ~ facet_col) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    strip.text = element_text(size = 18),         # Facet labels (e.g., column/row labels)
    axis.title = element_text(size = 18),         # Axis titles (x and y)
    axis.text = element_text(size = 14)
  ) +
  labs(
    title = "Parameters for choice-specific characteristics",
    x = "Time",
    y = "Beta_t"
  )


# Intercepts
indeces <- 10:12
# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
for (j in indeces) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = EP_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = EP_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j-9,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j-9,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = SMC_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = SMC_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(SMC_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = SMC_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(SMC_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+2],
    facet_row = j-9,
    facet_col = "SMC"
  )
  data_list[[length(data_list) + 1]] <- df_filter
}

plot_data <- bind_rows(data_list)

# Create the ggplot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), alpha = 0.5) +
  geom_line(aes(y = y2, color = color), alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(facet_row ~ facet_col) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    strip.text = element_text(size = 18),         # Facet labels (e.g., column/row labels)
    axis.title = element_text(size = 18),         # Axis titles (x and y)
    axis.text = element_text(size = 14)
  ) +
  labs(
    title = "Intercept parameters",
    x = "Time",
    y = "Beta_t"
    )


