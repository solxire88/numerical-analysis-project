# ============================================================
# Simple ML Project in R:
# Polynomial Linear Regression (4 coefficients) trained with
#   1) Jacobi (iterative)
#   2) Gauss–Seidel (iterative)
#   3) Relaxation / SOR (Successive Over-Relaxation)
#
# Model (4 coefficients):
#   y = w0 + w1*x + w2*x^2 + w3*x^3 + noise
#
# We train by solving the Normal Equation:
#   (X^T X) w = X^T y
# This is a linear system: A w = b
# ============================================================

set.seed(42)

# -----------------------------
# 1) Create sample data
# -----------------------------
n <- 140
x <- runif(n, -2, 2)

# True coefficients (the 4 unknowns we want to estimate)
w_true <- c(2.0, -1.5, 0.7, 0.2)  # (w0, w1, w2, w3)

# Design matrix X (n x 4): [1, x, x^2, x^3]
X <- cbind(1, x, x^2, x^3)

# Target y with noise
noise <- rnorm(n, mean = 0, sd = 0.4)
y <- as.numeric(X %*% w_true + noise)

# -----------------------------
# 2) Build the linear system A w = b
# -----------------------------
A <- t(X) %*% X
b <- t(X) %*% y

# IMPORTANT for iterative solvers:
# - Iterative methods converge more reliably if A is "nice" (e.g., SPD / diagonally dominant).
# - A = X^T X is symmetric and usually PSD. Adding a tiny ridge makes it SPD.
ridge <- 1e-8
A <- A + ridge * diag(nrow(A))

# ============================================================
# 3) ITERATIVE SOLVERS
#    All return a list: w, iters, residual_norm
# ============================================================

# Helper: residual norm ||A w - b||
residual_norm <- function(A, w, b) sqrt(sum((A %*% w - b)^2))

# -----------------------------
# 3.1 Jacobi Method
#   w^(k+1)_i = (b_i - sum_{j!=i} A_ij * w^(k)_j) / A_ii
# -----------------------------
jacobi_solve <- function(A, b, w0 = NULL, tol = 1e-10, max_iter = 5000) {
  A <- as.matrix(A); b <- as.numeric(b)
  n <- nrow(A)
  if (is.null(w0)) w0 <- rep(0, n)

  w_old <- w0
  D <- diag(A)                     # diagonal of A
  if (any(abs(D) < 1e-14)) stop("Jacobi: zero diagonal entry.")

  for (k in 1:max_iter) {
    w_new <- w_old
    for (i in 1:n) {
      # sum_{j!=i} A[i,j] * w_old[j]
      s <- sum(A[i, ] * w_old) - A[i, i] * w_old[i]
      w_new[i] <- (b[i] - s) / A[i, i]
    }
    if (sqrt(sum((w_new - w_old)^2)) < tol) {
      return(list(w = w_new, iters = k, residual = residual_norm(A, w_new, b)))
    }
    w_old <- w_new
  }
  list(w = w_old, iters = max_iter, residual = residual_norm(A, w_old, b))
}

# -----------------------------
# 3.2 Gauss–Seidel Method
#   Uses updated values immediately:
#   w_i <- (b_i - sum_{j<i} A_ij*w_j(new) - sum_{j>i} A_ij*w_j(old)) / A_ii
# -----------------------------
gauss_seidel_solve <- function(A, b, w0 = NULL, tol = 1e-10, max_iter = 5000) {
  A <- as.matrix(A); b <- as.numeric(b)
  n <- nrow(A)
  if (is.null(w0)) w0 <- rep(0, n)

  w <- w0
  if (any(abs(diag(A)) < 1e-14)) stop("Gauss-Seidel: zero diagonal entry.")

  for (k in 1:max_iter) {
    w_old <- w
    for (i in 1:n) {
      s1 <- if (i == 1) 0 else sum(A[i, 1:(i - 1)] * w[1:(i - 1)])     # new values
      s2 <- if (i == n) 0 else sum(A[i, (i + 1):n] * w_old[(i + 1):n]) # old values
      w[i] <- (b[i] - s1 - s2) / A[i, i]
    }
    if (sqrt(sum((w - w_old)^2)) < tol) {
      return(list(w = w, iters = k, residual = residual_norm(A, w, b)))
    }
  }
  list(w = w, iters = max_iter, residual = residual_norm(A, w, b))
}

# -----------------------------
# 3.3 Relaxation Method (SOR)
#   w_i <- (1-omega)*w_i(old) + omega*(GaussSeidel update)
#
#   omega = 1   -> Gauss-Seidel
#   1 < omega < 2 often speeds up (over-relaxation)
#   0 < omega < 1 under-relaxation (more stable sometimes)
# -----------------------------
sor_solve <- function(A, b, omega = 1.2, w0 = NULL, tol = 1e-10, max_iter = 5000) {
  A <- as.matrix(A); b <- as.numeric(b)
  n <- nrow(A)
  if (is.null(w0)) w0 <- rep(0, n)

  w <- w0
  if (any(abs(diag(A)) < 1e-14)) stop("SOR: zero diagonal entry.")
  if (omega <= 0 || omega >= 2) stop("SOR: choose omega in (0, 2).")

  for (k in 1:max_iter) {
    w_old <- w
    for (i in 1:n) {
      s1 <- if (i == 1) 0 else sum(A[i, 1:(i - 1)] * w[1:(i - 1)])     # new
      s2 <- if (i == n) 0 else sum(A[i, (i + 1):n] * w_old[(i + 1):n]) # old
      gs_update <- (b[i] - s1 - s2) / A[i, i]
      w[i] <- (1 - omega) * w_old[i] + omega * gs_update
    }
    if (sqrt(sum((w - w_old)^2)) < tol) {
      return(list(w = w, iters = k, residual = residual_norm(A, w, b)))
    }
  }
  list(w = w, iters = max_iter, residual = residual_norm(A, w, b))
}

# ============================================================
# 4) Solve for coefficients using the 3 methods
# ============================================================

w0_init <- rep(0, 4)

sol_jacobi <- jacobi_solve(A, b, w0 = w0_init)
sol_gs     <- gauss_seidel_solve(A, b, w0 = w0_init)
sol_sor    <- sor_solve(A, b, omega = 1.2, w0 = w0_init)

w_jacobi <- sol_jacobi$w
w_gs     <- sol_gs$w
w_sor    <- sol_sor$w

# -----------------------------
# 5) Display coefficients
# -----------------------------
terms <- c("w0 (intercept)", "w1 (x)", "w2 (x^2)", "w3 (x^3)")
coef_table <- data.frame(
  term  = terms,
  true  = w_true,
  jacobi = w_jacobi,
  gauss_seidel = w_gs,
  sor = w_sor
)
print(coef_table)

cat("\nIterations & residual norms:\n")
cat("Jacobi       :", sol_jacobi$iters, "iters | residual =", sol_jacobi$residual, "\n")
cat("Gauss-Seidel :", sol_gs$iters,     "iters | residual =", sol_gs$residual, "\n")
cat("SOR (omega=1.2):", sol_sor$iters,  "iters | residual =", sol_sor$residual, "\n")

# ============================================================
# 6) Plot the regression (data + fitted curves)
# ============================================================

# Smooth grid for curves
xg <- seq(min(x), max(x), length.out = 300)
Xg <- cbind(1, xg, xg^2, xg^3)

yg_j <- as.numeric(Xg %*% w_jacobi)
yg_g <- as.numeric(Xg %*% w_gs)
yg_s <- as.numeric(Xg %*% w_sor)

plot(x, y,
     pch = 16,
     main = "Polynomial Regression (4 coeffs) via Jacobi / Gauss-Seidel / SOR",
     xlab = "x", ylab = "y")

lines(xg, yg_j, lwd = 2, lty = 1)
lines(xg, yg_g, lwd = 2, lty = 2)
lines(xg, yg_s, lwd = 2, lty = 3)

legend("topleft",
       legend = c("Jacobi", "Gauss-Seidel", "SOR (Relaxation)"),
       lty = c(1, 2, 3), lwd = 2, bty = "n")
