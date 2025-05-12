# Load necessary package
library(Matrix)  # for matrix utilities (optional)

# Define fixed parameters
m <- 2  # any positive scalar
V <- matrix(c(2, 0.5, 0.3,
              0.5, 1.5, 0.4,
              0.3, 0.4, 1), nrow = 3, byrow = TRUE)
I <- diag(3)

# Function that returns A as a function of c
A_func <- function(c) {
  (1/m) * (V - 2 * c * I)
}

# Function to compute trace of inverse of A(c)
trace_inv_A <- function(c) {
  A <- A_func(c)
  
  if (det(A) <= 1e-10) {
    return(Inf)  # Avoid singular or near-singular matrices
  }
  
  A_inv <- solve(A)
  sum(diag(A_inv))
}

# Root-finding function: solve trace_inv_A(c) = 5
target_trace <- 5
root_function <- function(c) {
  trace_inv_A(c) - target_trace
}

# Use uniroot to find c in a reasonable interval
# You might need to adjust the interval depending on V and m
result <- uniroot(root_function, lower = -1, upper = 2)

# Show the result
cat("Value of c such that trace(A^-1) = 5:", result$root, "\n")
