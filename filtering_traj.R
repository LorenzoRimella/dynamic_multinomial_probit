library(dplyr)

library(ggplot2)

# T = 10, lambda_m = 10
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_10/filtering_trajectory_p_9_mt_10_L_3_Time_10_data_multinomial_probit_exact_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_10/filtering_trajectory_p_9_mt_10_L_3_Time_10_data_multinomial_probit_SMC_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_10/filtering_trajectory_p_9_mt_10_L_3_Time_10_data_multinomial_probit_ADF_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_10/filtering_trajectory_p_9_mt_10_L_3_Time_10_data_multinomial_probit_EP_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_10/filtering_trajectory_p_9_mt_10_L_3_Time_10_data_multinomial_probit_single.Rdata")

Time <- length(data_multinomial_probit$Y)

color_list <- c("blue", "black", "red", "magenta", "brown")

# individual-specific characteristics
indeces <- c(1, 2, 3, 4, 5)
# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
i=1
for (j in indeces) {
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = exact_filter_sampler_output$a_history[, j],
    y_upper = exact_filter_sampler_output$a_history[, j] + 1.96 * sqrt(exact_filter_sampler_output$P_history[, j, j]),
    y_lower = exact_filter_sampler_output$a_history[, j] - 1.96 * sqrt(exact_filter_sampler_output$P_history[, j, j]),
    color = color_list[1],
    facet_row = i,
    facet_col = "i.i.d."
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = EP_filter_sampler_output$a_history[, j],
    y_upper = EP_filter_sampler_output$a_history[, j] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, j, j]),
    y_lower = EP_filter_sampler_output$a_history[, j] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, j, j]),
    color = color_list[2],
    facet_row = i,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = ADF_filter_sampler_output$a_history[, j],
    y_upper = ADF_filter_sampler_output$a_history[, j] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, j, j]),
    y_lower = ADF_filter_sampler_output$a_history[, j] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, j, j]),
    color = color_list[3],
    facet_row = i,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = SMC_filter_sampler_output$a_history[, j],
    y_upper = SMC_filter_sampler_output$a_history[, j] + 1.96 * sqrt(SMC_filter_sampler_output$P_history[, j, j]),
    y_lower = SMC_filter_sampler_output$a_history[, j] - 1.96 * sqrt(SMC_filter_sampler_output$P_history[, j, j]),
    color = color_list[4],
    facet_row = i,
    facet_col = "SMC"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  i=i+1
}

plot_data <- bind_rows(data_list)

# specify order of the methods
plot_data$facet_col = factor(plot_data$facet_col, levels = c("i.i.d.", "ADF", "EP", "SMC"))

# Create the ggplot
gg1 = ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), linetype = "dashed", alpha = 0.8) +
  geom_line(aes(y = y2, color = color), linetype = "solid", alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(
    facet_row ~ facet_col,
    labeller = labeller(
      facet_row = setNames(
        lapply(indeces, \(k) bquote(beta[t]^.(k))),
        indeces
      )
    )
  )+
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14),         # Facet labels (e.g., column/row labels)
    axis.title = element_text(size = 14),         # Axis titles (x and y)
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "",
    x = "Time",
    y = expression(beta[t])
  )

ggsave("figures/traj_T10_m10.pdf", plot = gg1, width = 10, height = 6, units = "in")
ggsave("figures/traj_T10_m10.png", plot = gg1, width = 10, height = 6, dpi = 300, units = "in")


# T = 50, lambda_m = 10
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_50/filtering_trajectory_p_9_mt_10_L_3_Time_50_data_multinomial_probit_SMC_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_50/filtering_trajectory_p_9_mt_10_L_3_Time_50_data_multinomial_probit_ADF_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_50/filtering_trajectory_p_9_mt_10_L_3_Time_50_data_multinomial_probit_EP_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_10_L_3_Time_50/filtering_trajectory_p_9_mt_10_L_3_Time_50_data_multinomial_probit_single.Rdata")

Time <- length(data_multinomial_probit$Y)

# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
i=1
for (j in indeces) {
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = EP_filter_sampler_output$a_history[, j],
    y_upper = EP_filter_sampler_output$a_history[, j] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, j, j]),
    y_lower = EP_filter_sampler_output$a_history[, j] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, j, j]),
    color = color_list[2],
    facet_row = i,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = ADF_filter_sampler_output$a_history[, j],
    y_upper = ADF_filter_sampler_output$a_history[, j] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, j, j]),
    y_lower = ADF_filter_sampler_output$a_history[, j] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, j, j]),
    color = color_list[3],
    facet_row = i,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = SMC_filter_sampler_output$a_history[, j],
    y_upper = SMC_filter_sampler_output$a_history[, j] + 1.96 * sqrt(SMC_filter_sampler_output$P_history[, j, j]),
    y_lower = SMC_filter_sampler_output$a_history[, j] - 1.96 * sqrt(SMC_filter_sampler_output$P_history[, j, j]),
    color = color_list[4],
    facet_row = i,
    facet_col = "SMC"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  i=i+1
}

plot_data <- bind_rows(data_list)

# specify order of the methods
plot_data$facet_col = factor(plot_data$facet_col, levels = c("i.i.d.", "ADF", "EP", "SMC"))

# Create the ggplot
gg2 = ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), linetype = "dashed", alpha = 0.8) +
  geom_line(aes(y = y2, color = color), linetype = "solid", alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(
    facet_row ~ facet_col,
    labeller = labeller(
      facet_row = setNames(
        lapply(indeces, \(k) bquote(beta[t]^.(k))),
        indeces
      )
    )
  )+
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14),         # Facet labels (e.g., column/row labels)
    axis.title = element_text(size = 14),         # Axis titles (x and y)
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "",
    x = "Time",
    y = expression(beta[t])
  )


ggsave("figures/traj_T50_m10.pdf", plot = gg2, width = 10, height = 6, units = "in")
ggsave("figures/traj_T50_m10.png", plot = gg2, width = 10, height = 6, dpi = 300, units = "in")

# T = 50, lambda_m = 50
load("data/filtering_trajectory_p_9_mt_50_L_3_Time_50/filtering_trajectory_p_9_mt_50_L_3_Time_50_data_multinomial_probit_ADF_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_50_L_3_Time_50/filtering_trajectory_p_9_mt_50_L_3_Time_50_data_multinomial_probit_EP_filter.Rdata")
load("data/filtering_trajectory_p_9_mt_50_L_3_Time_50/filtering_trajectory_p_9_mt_50_L_3_Time_50_data_multinomial_probit_single.Rdata")

# PLOT FILTERING
# Create a data frame for plotting

Time = 50
data_list <- list()
i=1
for (j in indeces) {
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = EP_filter_sampler_output$a_history[, j],
    y_upper = EP_filter_sampler_output$a_history[, j] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, j, j]),
    y_lower = EP_filter_sampler_output$a_history[, j] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, j, j]),
    color = color_list[2],
    facet_row = i,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, j],
    y2 = ADF_filter_sampler_output$a_history[, j],
    y_upper = ADF_filter_sampler_output$a_history[, j] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, j, j]),
    y_lower = ADF_filter_sampler_output$a_history[, j] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, j, j]),
    color = color_list[3],
    facet_row = i,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  i=i+1
}

plot_data <- bind_rows(data_list)

# specify order of the methods
plot_data$facet_col = factor(plot_data$facet_col, levels = c("i.i.d.", "ADF", "EP", "SMC"))

# Create the ggplot
gg3 = ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), linetype = "dashed", alpha = 0.8) +
  geom_line(aes(y = y2, color = color), linetype = "solid", alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(
    facet_row ~ facet_col,
    labeller = labeller(
      facet_row = setNames(
        lapply(indeces, \(k) bquote(beta[t]^.(k))),
        indeces
      )
    )
  )+
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14),         # Facet labels (e.g., column/row labels)
    axis.title = element_text(size = 14),         # Axis titles (x and y)
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "",
    x = "Time",
    y = expression(beta[t])
  )


ggsave("figures/traj_T50_m50.pdf", plot = gg3, width = 10, height = 6, units = "in")
ggsave("figures/traj_T50_m50.png", plot = gg3, width = 10, height = 6, dpi = 300, units = "in")
