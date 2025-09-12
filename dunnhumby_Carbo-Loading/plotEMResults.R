load("CarboLoading_EM_output_T30_200.RData")

plot(EM_output$log_likelihood_history, type = "l")

str(EM_output$Sigma_history)


# A: 3D array with dims c(I, 5, 5) so that A[i,,] is 5x5
# Set include_diag = TRUE if you also want the diagonal
plot_lower_tri_trajectories <- function(A, include_diag = FALSE) {
  stopifnot(length(dim(A)) == 3, dim(A)[2] == 5, dim(A)[3] == 5)
  I <- dim(A)[1]
  
  # Lower-triangular index (toggle diag with include_diag)
  lt <- lower.tri(matrix(0, 5, 5), diag = include_diag)
  rc <- which(lt, arr.ind = TRUE)
  
  # Optional: use dimnames if available for nicer labels
  rnames <- dimnames(A)[[2]]
  cnames <- dimnames(A)[[3]]
  if (is.null(rnames)) rnames <- as.character(1:5)
  if (is.null(cnames)) cnames <- as.character(1:5)
  
  pair_lab <- paste0(rnames[rc[,1]], " vs ", cnames[rc[,2]])
  
  # Collect values for each (row, col) across i
  # vals: I x K matrix (K = number of lower-tri pairs)
  vals <- sapply(seq_len(nrow(rc)), function(k) A[, rc[k,1], rc[k,2]])
  
  df <- data.frame(
    i = rep(seq_len(I), times = ncol(vals)),
    pair = rep(pair_lab, each = I),
    value = as.vector(vals)
  )
  
  # Plot
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2")
  }
  library(ggplot2)
  ggplot(df, aes(x = i, y = value, color = pair)) +
    geom_line() +
    labs(title = "Lower-triangular trajectories", x = "i", y = "A[i, row, col]", color = "Pair") +
    theme_minimal()
}

plot_lower_tri_trajectories(EM_output$Sigma_history, include_diag = T)
