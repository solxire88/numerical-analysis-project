# ============================================================
# Simple Linear Regression (4 coefficients) solved 3 ways in R
#   1) Méthode de Gauss (Gaussian elimination)
#   2) Décomposition LU (Doolittle + pivot partiel)
#   3) Décomposition de Cholesky (A = L L^T)
#
# Model (4 coefficients):
#   y = w0 + w1*x + w2*x^2 + w3*x^3 + noise
#
# Training uses the Normal Equation:
#   (X^T X) w = X^T y
# where X = [1, x, x^2, x^3]
# ============================================================

set.seed(42)

# -----------------------------
# 1) Create sample data
# -----------------------------
n <- 120
x <- runif(n, -2, 2)

# True coefficients (the 4 coefficients we want to recover)
w_true <- c(2.0, -1.5, 0.7, 0.2)  # (w0, w1, w2, w3)

# Build design matrix X (n x 4): [1, x, x^2, x^3]
X <- cbind(1, x, x^2, x^3)

# Generate target y with noise
noise <- rnorm(n, mean = 0, sd = 0.4)
y <- as.numeric(X %*% w_true + noise)

# -----------------------------
# 2) Normal equation matrices
#    A w = b
# -----------------------------
A <- t(X) %*% X
b <- t(X) %*% y

# Optional: tiny ridge term for stability (especially for Cholesky)
# You can set ridge <- 0 to keep it "pure".
ridge <- 1e-10
A_ridge <- A + ridge * diag(nrow(A))

# ============================================================
# 3) SOLVERS
# ============================================================

# -----------------------------
# 3.1 Gaussian Elimination (Gauss) with partial pivoting
# -----------------------------
gauss_solve <- function(A, b, tol = 1e-12) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  Aug <- cbind(A, b)  # augmented matrix [A | b]
  
  # Forward elimination
  for (k in 1:(n - 1)) {
    # Partial pivot: choose row with largest absolute pivot
    pivot_row <- k - 1 + which.max(abs(Aug[k:n, k]))
    if (abs(Aug[pivot_row, k]) < tol) stop("Singular / nearly singular matrix in Gauss.")
    
    # Swap rows if needed
    if (pivot_row != k) {
      tmp <- Aug[k, ]; Aug[k, ] <- Aug[pivot_row, ]; Aug[pivot_row, ] <- tmp
    }
    
    # Eliminate entries below the pivot
    for (i in (k + 1):n) {
      factor <- Aug[i, k] / Aug[k, k]
      Aug[i, k:(n + 1)] <- Aug[i, k:(n + 1)] - factor * Aug[k, k:(n + 1)]
    }
  }
  
  # Back substitution
  w <- numeric(n)
  w[n] <- Aug[n, n + 1] / Aug[n, n]
  for (i in (n - 1):1) {
    w[i] <- (Aug[i, n + 1] - sum(Aug[i, (i + 1):n] * w[(i + 1):n])) / Aug[i, i]
  }
  w
}

# -----------------------------
# 3.2 LU Decomposition (Doolittle) + partial pivoting
# -----------------------------
lu_decomp <- function(A, tol = 1e-12) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  U <- A
  L <- diag(1, n)
  P <- diag(1, n)
  
  for (k in 1:(n - 1)) {
    # pivot row for column k
    pivot_row <- k - 1 + which.max(abs(U[k:n, k]))
    if (abs(U[pivot_row, k]) < tol) stop("Singular / nearly singular matrix in LU.")
    
    # swap rows in U and P; also swap in L up to column k-1
    if (pivot_row != k) {
      tmp <- U[k, ]; U[k, ] <- U[pivot_row, ]; U[pivot_row, ] <- tmp
      tmp <- P[k, ]; P[k, ] <- P[pivot_row, ]; P[pivot_row, ] <- tmp
      if (k > 1) {
        tmp <- L[k, 1:(k - 1)]
        L[k, 1:(k - 1)] <- L[pivot_row, 1:(k - 1)]
        L[pivot_row, 1:(k - 1)] <- tmp
      }
    }
    
    # elimination to form L and U
    for (i in (k + 1):n) {
      L[i, k] <- U[i, k] / U[k, k]
      U[i, k:n] <- U[i, k:n] - L[i, k] * U[k, k:n]
    }
  }
  
  list(L = L, U = U, P = P)
}

forward_sub <- function(L, b, tol = 1e-12) {
  n <- nrow(L)
  y <- numeric(n)
  for (i in 1:n) {
    if (abs(L[i, i]) < tol) stop("Zero diagonal in forward substitution.")
    s <- if (i == 1) 0 else sum(L[i, 1:(i - 1)] * y[1:(i - 1)])
    y[i] <- (b[i] - s) / L[i, i]
  }
  y
}

back_sub <- function(U, y, tol = 1e-12) {
  n <- nrow(U)
  x <- numeric(n)
  for (i in n:1) {
    if (abs(U[i, i]) < tol) stop("Zero diagonal in back substitution.")
    s <- if (i == n) 0 else sum(U[i, (i + 1):n] * x[(i + 1):n])
    x[i] <- (y[i] - s) / U[i, i]
  }
  x
}

lu_solve <- function(A, b) {
  dec <- lu_decomp(A)
  # Solve P*A*w = P*b  =>  L*U*w = P*b
  Pb <- as.numeric(dec$P %*% b)
  y <- forward_sub(dec$L, Pb)
  w <- back_sub(dec$U, y)
  w
}

# -----------------------------
# 3.3 Cholesky Decomposition (A = L L^T)
#     Requires A to be symmetric positive definite (SPD).
#     Here A = X^T X is symmetric; ridge helps guarantee SPD.
# -----------------------------
cholesky_decomp <- function(A, tol = 1e-12) {
  A <- as.matrix(A)
  n <- nrow(A)
  L <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:i) {
      s <- if (j == 1) 0 else sum(L[i, 1:(j - 1)] * L[j, 1:(j - 1)])
      
      if (i == j) {
        val <- A[i, i] - s
        if (val <= tol) stop("Cholesky failed: A not positive definite.")
        L[i, j] <- sqrt(val)
      } else {
        L[i, j] <- (A[i, j] - s) / L[j, j]
      }
    }
  }
  L
}

cholesky_solve <- function(A, b) {
  L <- cholesky_decomp(A)
  # Solve L z = b, then L^T w = z
  z <- forward_sub(L, b)
  w <- back_sub(t(L), z)
  w
}

# ============================================================
# 4) Solve for coefficients with the 3 methods
# ============================================================

w_gauss <- gauss_solve(A, b)
w_lu    <- lu_solve(A, b)
w_chol  <- cholesky_solve(A_ridge, b)  # using ridge-stabilized A

# -----------------------------
# 5) Display coefficients
# -----------------------------
terms <- c("w0 (intercept)", "w1 (x)", "w2 (x^2)", "w3 (x^3)")
coef_table <- data.frame(
  term  = terms,
  true  = w_true,
  gauss = w_gauss,
  lu    = w_lu,
  chol  = w_chol
)
print(coef_table)

# ============================================================
# 6) Plot: data + fitted regression curve
# ============================================================

# Create a smooth x-grid to draw curves
xg <- seq(min(x), max(x), length.out = 300)
Xg <- cbind(1, xg, xg^2, xg^3)

yg_gauss <- as.numeric(Xg %*% w_gauss)
yg_lu    <- as.numeric(Xg %*% w_lu)
yg_chol  <- as.numeric(Xg %*% w_chol)

# Scatter plot of data
plot(x, y,
     pch = 16,
     main = "Polynomial Regression (4 coefficients) solved by Gauss / LU / Cholesky",
     xlab = "x",
     ylab = "y")

# Regression curves (they should almost overlap)
lines(xg, yg_gauss, lwd = 2, lty = 1)
lines(xg, yg_lu,    lwd = 2, lty = 2)
lines(xg, yg_chol,  lwd = 2, lty = 3)

legend("topleft",
       legend = c("Gauss", "LU", "Cholesky"),
       lty = c(1, 2, 3),
       lwd = 2,
       bty = "n")
