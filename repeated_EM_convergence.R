library(ggplot2)

# This script reads the saved EM output without rerunning the experiment.
times <- c(10L, 30L, 60L)
rolling_window <- 10L
tail_window <- 20L

result_file <- function(Time) {
  prefix <- sprintf(
    "Sigma_recovery_p_9_mt_15_L_3_Time_%d",
    Time
  )
  file.path(
    "data",
    prefix,
    paste0(prefix, "_log_likelihood_history.Rdata")
  )
}

rolling_mean <- function(x, window) {
  as.numeric(
    stats::filter(
      x,
      rep(1 / window, window),
      sides = 1
    )
  )
}

diagnostic_data <- lapply(times, function(Time) {
  input_file <- result_file(Time)
  if (!file.exists(input_file)) {
    stop("Saved log-likelihood history not found: ", input_file)
  }

  saved_output <- new.env(parent = emptyenv())
  load(input_file, envir = saved_output)
  loglik <- saved_output$log_likelihood_history

  if (!is.matrix(loglik) || any(!is.finite(loglik))) {
    stop("The saved log-likelihood history is not a complete numeric matrix: ",
         input_file)
  }

  iterations <- seq.int(0L, ncol(loglik) - 1L)
  smoothed <- t(apply(loglik, 1L, rolling_mean, window = rolling_window))

  # Different replications use different simulated responses. Their absolute
  # likelihood levels are therefore not comparable. Center each trajectory at
  # its own average smoothed likelihood over the final iterations.
  tail_columns <- seq.int(ncol(smoothed) - tail_window + 1L, ncol(smoothed))
  tail_reference <- rowMeans(smoothed[, tail_columns, drop = FALSE])
  centered <- sweep(smoothed, 1L, tail_reference, FUN = "-")

  quantiles <- apply(
    centered,
    2L,
    stats::quantile,
    probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    na.rm = TRUE,
    names = FALSE
  )

  data.frame(
    Time = sprintf("T = %d", Time),
    iteration = iterations,
    q05 = quantiles[1L, ],
    q25 = quantiles[2L, ],
    median = quantiles[3L, ],
    q75 = quantiles[4L, ],
    q95 = quantiles[5L, ],
    n_replications = nrow(loglik)
  )
})

diagnostic_data <- do.call(rbind, diagnostic_data)
diagnostic_data <- diagnostic_data[is.finite(diagnostic_data$median), ]
diagnostic_data$Time <- factor(
  diagnostic_data$Time,
  levels = sprintf("T = %d", times)
)

convergence_plot <- ggplot(
  diagnostic_data,
  aes(x = iteration, y = median)
) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_ribbon(
    aes(ymin = q05, ymax = q95),
    fill = "#7FA6CC",
    alpha = 0.22
  ) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75),
    fill = "#3975A8",
    alpha = 0.32
  ) +
  geom_line(color = "#174A73", linewidth = 0.75) +
  facet_wrap(~Time, nrow = 1L, scales = "free_y") +
  scale_x_continuous(
    breaks = c(10, 50, 100, 150, 200),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "EM iteration",
    y = "Centered smoothed log likelihood"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    panel.spacing = grid::unit(0.8, "lines"),
    strip.text = element_text(face = "bold", color = "black"),
    plot.margin = margin(5, 7, 5, 5)
  )

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave(
  "figures/em_repeated_convergence.pdf",
  plot = convergence_plot,
  width = 7.2,
  height = 2.7,
  units = "in",
  bg = "white"
)
ggsave(
  "figures/em_repeated_convergence.png",
  plot = convergence_plot,
  width = 7.2,
  height = 2.7,
  units = "in",
  dpi = 320,
  bg = "white"
)

# A scale-free numerical check: compare the change between the last two
# 25-iteration blocks with the total improvement from the first 25 iterations.
tail_diagnostic <- lapply(times, function(Time) {
  saved_output <- new.env(parent = emptyenv())
  load(result_file(Time), envir = saved_output)
  loglik <- saved_output$log_likelihood_history

  n_iterations_saved <- ncol(loglik)
  early_mean <- rowMeans(loglik[, 1L:25L, drop = FALSE])
  penultimate_columns <- seq.int(
    n_iterations_saved - 49L,
    n_iterations_saved - 25L
  )
  final_columns <- seq.int(
    n_iterations_saved - 24L,
    n_iterations_saved
  )
  penultimate_mean <- rowMeans(
    loglik[, penultimate_columns, drop = FALSE]
  )
  final_mean <- rowMeans(loglik[, final_columns, drop = FALSE])
  relative_tail_change <- abs(final_mean - penultimate_mean) /
    pmax(abs(final_mean - early_mean), .Machine$double.eps)

  data.frame(
    Time = Time,
    replications = nrow(loglik),
    median_relative_tail_change = stats::median(relative_tail_change),
    q90_relative_tail_change = unname(
      stats::quantile(relative_tail_change, 0.90)
    )
  )
})

tail_diagnostic <- do.call(rbind, tail_diagnostic)
print(tail_diagnostic, row.names = FALSE, digits = 3)

