library(dplyr)
library(xtable)

summarise_filtering_results <- function(file, reference_type) {
  saved_output <- new.env(parent = emptyenv())
  load(file, envir = saved_output)
  results <- saved_output$data_frame_output

  if (!"replicate" %in% names(results)) {
    stop("The results do not contain the replicate identifier: ", file)
  }

  reference_loglik <- results %>%
    filter(type == reference_type) %>%
    select(replicate, reference_log_likelihood = log_likelihood)

  results <- results %>%
    left_join(reference_loglik, by = "replicate") %>%
    transmute(
      replicate,
      type,
      filter_time = filter_running_time,
      delta_log_likelihood = log_likelihood - reference_log_likelihood,
      mean_marginal_W2 = filter_distance_from_exact
    )

  numeric_summary <- results %>%
    group_by(type) %>%
    summarise(
      across(
        c(filter_time, delta_log_likelihood, mean_marginal_W2),
        list(mean = mean, sd = sd),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  formatted_summary <- numeric_summary %>%
    transmute(
      type,
      filter_time = sprintf("%.3f ± %.3f", filter_time_mean, filter_time_sd),
      delta_log_likelihood = sprintf(
        "%.3f ± %.3f",
        delta_log_likelihood_mean,
        delta_log_likelihood_sd
      ),
      mean_marginal_W2 = sprintf(
        "%.3f ± %.3f",
        mean_marginal_W2_mean,
        mean_marginal_W2_sd
      )
    )

  print(formatted_summary, row.names = FALSE)
  print(xtable(formatted_summary), include.rownames = FALSE)
  invisible(numeric_summary)
}

summary_T10_m10 <- summarise_filtering_results(
  paste0(
    "data/full_parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_10/",
    "full_parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_10_data.Rdata"
  ),
  reference_type = "exact"
)

summary_T50_m10 <- summarise_filtering_results(
  paste0(
    "data/full_parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_50/",
    "full_parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_50_data.Rdata"
  ),
  reference_type = "ADF"
)

summary_T50_m50 <- summarise_filtering_results(
  paste0(
    "data/full_parallel_filt_smoot_lik_p_9_mt_50_L_3_Time_50/",
    "full_parallel_filt_smoot_lik_p_9_mt_50_L_3_Time_50_data.Rdata"
  ),
  reference_type = "ADF"
)
