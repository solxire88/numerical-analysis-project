############################################################
# Simple Imputation (Interpolation) in R
# - Time series example: Temperature over 7 days
# - Create a "hole" (missing value) in the middle
# - Impute the missing value using:
#     1) Direct (naive) interpolation (linear between neighbors)
#     2) Lagrange interpolation polynomial
#     3) Newton interpolation (divided differences)
# - Compare Imputed vs Real (hidden) value
# - Plot the series + interpolation curves filling the gap
#
# Based on standard polynomial interpolation ideas:
# - Direct/naive method (Vandermonde / simple approach)
# - Lagrange polynomial
# - Newton divided differences
# (see your course chapter on interpolation) :contentReference[oaicite:0]{index=0}
############################################################

set.seed(42)

#############################
# 1) Create sample time series
#############################
days <- 1:7

# Create a smooth-ish temperature pattern with small noise
temp_true <- 20 + 2*sin(days/2) + rnorm(length(days), 0, 0.15)

# Create a missing value ("hole") in the middle: day 4
missing_day <- 4
real_hidden <- temp_true[missing_day]

temp_obs <- temp_true
temp_obs[missing_day] <- NA

cat("Real hidden value at day", missing_day, "=", real_hidden, "\n")

###############################################
# Helper: pick points around the missing day
# We'll use 4 surrounding points -> degree 3 poly:
# days = 2,3,5,6 (avoid using the missing point)
###############################################
idx_used <- c(2, 3, 5, 6)
x_known <- days[idx_used]
y_known <- temp_obs[idx_used]
x_star <- missing_day  # where we want to impute

#############################
# 2) Method 1: Direct (naive) - linear interpolation
# Use immediate neighbors day 3 and day 5
#############################
x0 <- 3; x1 <- 5
y0 <- temp_obs[x0]; y1 <- temp_obs[x1]

# Linear interpolation formula:
# y(x) = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
impute_direct <- y0 + (y1 - y0) * (x_star - x0) / (x1 - x0)

#############################
# 3) Method 2: Lagrange interpolation
# P(x) = Σ y_i * L_i(x)
# where L_i(x)= Π_{j!=i} (x - x_j)/(x_i - x_j)
#############################
lagrange_eval <- function(x, y, x_star) {
  n <- length(x)
  p <- 0
  for (i in 1:n) {
    Li <- 1
    for (j in 1:n) {
      if (j != i) {
        Li <- Li * (x_star - x[j]) / (x[i] - x[j])
      }
    }
    p <- p + y[i] * Li
  }
  p
}

impute_lagrange <- lagrange_eval(x_known, y_known, x_star)

#############################
# 4) Method 3: Newton interpolation (divided differences)
# Build coefficients a_k = f[x0,...,xk]
# P(x)= a0 + a1(x-x0) + a2(x-x0)(x-x1)+...
#############################
newton_divided_diff <- function(x, y) {
  n <- length(x)
  dd <- matrix(0, n, n)
  dd[, 1] <- y

  for (j in 2:n) {
    for (i in 1:(n - j + 1)) {
      dd[i, j] <- (dd[i + 1, j - 1] - dd[i, j - 1]) / (x[i + j - 1] - x[i])
    }
  }
  dd[1, ]  # coefficients a0..a_{n-1}
}

newton_eval <- function(x, coef, x_star) {
  n <- length(coef)
  p <- coef[1]
  prod_term <- 1
  for (k in 2:n) {
    prod_term <- prod_term * (x_star - x[k - 1])
    p <- p + coef[k] * prod_term
  }
  p
}

coef_newton <- newton_divided_diff(x_known, y_known)
impute_newton <- newton_eval(x_known, coef_newton, x_star)

#############################
# 5) Display results
#############################
results <- data.frame(
  Method = c("Direct (Linear)", "Lagrange (deg 3)", "Newton (deg 3)"),
  Imputed = c(impute_direct, impute_lagrange, impute_newton),
  Real = rep(real_hidden, 3),
  Abs_Error = abs(c(impute_direct, impute_lagrange, impute_newton) - real_hidden)
)

cat("\n=== Imputation Results ===\n")
print(results)

#############################
# 6) Plot: time series + filled gap curves
#############################

# Create a fine grid to draw interpolation curves between day 2 and day 6
grid_x <- seq(2, 6, length.out = 200)

# Direct (piecewise linear between day 3 and day 5) for plotting
direct_curve <- y0 + (y1 - y0) * (grid_x - x0) / (x1 - x0)

# Lagrange curve using our 4 points
lagrange_curve <- sapply(grid_x, function(xx) lagrange_eval(x_known, y_known, xx))

# Newton curve using same 4 points
newton_curve <- sapply(grid_x, function(xx) newton_eval(x_known, coef_newton, xx))

# Plot observed series (with hole)
plot(days, temp_obs, type = "b", pch = 16,
     xlab = "Day", ylab = "Temperature",
     main = "Imputation by Interpolation (Direct vs Lagrange vs Newton)",
     ylim = range(c(temp_true, lagrange_curve, newton_curve), na.rm = TRUE))

# Mark the real hidden point (what we removed)
points(missing_day, real_hidden, pch = 8, cex = 1.6)  # star-like marker

# Mark the imputed values
points(missing_day, impute_direct,   pch = 1, cex = 1.4)
points(missing_day, impute_lagrange, pch = 2, cex = 1.4)
points(missing_day, impute_newton,   pch = 0, cex = 1.4)

# Draw the interpolation curves over the gap zone
lines(grid_x, direct_curve, lty = 2, lwd = 2)
lines(grid_x, lagrange_curve, lty = 1, lwd = 2)
lines(grid_x, newton_curve, lty = 3, lwd = 2)

legend("topleft",
       legend = c("Observed (with hole)", "Real hidden value", "Direct impute",
                  "Lagrange curve", "Newton curve"),
       pch = c(16, 8, 1, NA, NA),
       lty = c(NA, NA, NA, 1, 3),
       lwd = c(NA, NA, NA, 2, 2),
       bty = "n")
