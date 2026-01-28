############################################################
# SIMPLE PCA (Dimensionality Reduction) in R
# - Create highly correlated data (compresses well)
# - Compute PC1 using multiple eigen-iterative methods:
#     1) Power Iteration (puissance itérée)
#     2) Deflation (déflation)  [uses Power for PC1 then deflates for PC2 demo]
#     3) Jacobi rotations (Jacobi + rotations de Givens style for symmetric matrices)
#     4) QR iteration with Gram–Schmidt QR (fixed to avoid R indexing bug)
#     5) QR iteration with Householder QR
# - Reconstruct ("interpolate") data using only PC1 (rank-1 approximation)
# - Display Real vs Interpolated values
# - Plot PCA space + interpolated (reconstructed) points
############################################################

set.seed(42)

#############################
# 1) Generate correlated data
#############################
n <- 120
z <- rnorm(n)  # latent factor
X <- cbind(
  x1 = z + 0.20*rnorm(n),
  x2 = 0.9*z + 0.20*rnorm(n),
  x3 = -0.7*z + 0.20*rnorm(n)
)

# Standardize (center + scale) for PCA
Xs <- scale(X, center = TRUE, scale = TRUE)

# Covariance matrix (symmetric)
S <- cov(Xs)

# Small ridge for numerical stability (helps iterative methods & SPD)
S <- S + 1e-10 * diag(nrow(S))

#############################
# Helpers
#############################
norm2 <- function(v) sqrt(sum(v^2))

align_sign <- function(v, ref) {
  # eigenvectors can flip sign; align for fair comparison
  if (sum(v * ref) < 0) -v else v
}

rmse <- function(A, B) sqrt(mean((A - B)^2))

# Rank-1 reconstruction using PC1 vector v1:
# scores t = X v1 ; reconstructed Xhat = t v1^T
reconstruct_pc1 <- function(Xs, v1) {
  t1 <- Xs %*% v1
  Xhat <- t1 %*% t(v1)
  list(scores = as.vector(t1), Xhat = Xhat)
}

#############################
# Reference PCA (built-in eigen)
#############################
eig_ref <- eigen(S)
v1_true <- eig_ref$vectors[, 1]
v2_true <- eig_ref$vectors[, 2]
lambda1_true <- eig_ref$values[1]

############################################################
# 2) METHOD 1: Power Iteration (Puissance itérée)
############################################################
power_iteration <- function(A, max_iter = 2000, tol = 1e-12) {
  n <- nrow(A)
  v <- rnorm(n)
  v <- v / norm2(v)

  lam_old <- 0
  for (k in 1:max_iter) {
    w <- A %*% v
    v_new <- w / norm2(w)
    lam_new <- as.numeric(t(v_new) %*% A %*% v_new)  # Rayleigh quotient

    if (abs(lam_new - lam_old) < tol) break
    v <- v_new
    lam_old <- lam_new
  }
  list(value = lam_old, vector = v)
}

############################################################
# 3) METHOD 2: Deflation (Déflation)
############################################################
deflation_power <- function(A, max_iter = 2000, tol = 1e-12) {
  e1 <- power_iteration(A, max_iter, tol)
  v1 <- e1$vector
  lam1 <- e1$value

  # Deflate matrix (remove first component)
  A1 <- A - lam1 * (v1 %*% t(v1)) / as.numeric(t(v1) %*% v1)

  e2 <- power_iteration(A1, max_iter, tol)
  list(lam1 = lam1, v1 = v1, lam2 = e2$value, v2 = e2$vector)
}

############################################################
# 4) METHOD 3: Jacobi Eigen (rotations) for symmetric matrices
############################################################
jacobi_eigen <- function(A, max_iter = 5000, tol = 1e-12) {
  n <- nrow(A)
  V <- diag(n)

  # helper: find indices (p,q) of largest off-diagonal element
  max_offdiag <- function(M) {
    M2 <- abs(M)
    diag(M2) <- 0
    idx <- which(M2 == max(M2), arr.ind = TRUE)[1, ]
    list(p = idx[1], q = idx[2], val = M2[idx[1], idx[2]])
  }

  for (k in 1:max_iter) {
    m <- max_offdiag(A)
    if (m$val < tol) break
    p <- m$p; q <- m$q

    app <- A[p,p]; aqq <- A[q,q]; apq <- A[p,q]
    theta <- 0.5 * atan2(2 * apq, (app - aqq))
    c <- cos(theta); s <- sin(theta)

    # rotate rows/cols p and q (keep symmetry)
    for (i in 1:n) {
      if (i != p && i != q) {
        aip <- A[i,p]; aiq <- A[i,q]
        A[i,p] <- c*aip - s*aiq
        A[p,i] <- A[i,p]
        A[i,q] <- s*aip + c*aiq
        A[q,i] <- A[i,q]
      }
    }

    # update diagonal + zero out A[p,q]
    A[p,p] <- c^2*app - 2*s*c*apq + s^2*aqq
    A[q,q] <- s^2*app + 2*s*c*apq + c^2*aqq
    A[p,q] <- 0
    A[q,p] <- 0

    # update eigenvectors
    for (i in 1:n) {
      vip <- V[i,p]; viq <- V[i,q]
      V[i,p] <- c*vip - s*viq
      V[i,q] <- s*vip + c*viq
    }
  }

  list(values = diag(A), vectors = V)
}

############################################################
# 5) METHOD 4: QR Iteration with Gram–Schmidt QR (fixed)
############################################################
qr_gs <- function(A, tol = 1e-14) {
  A <- as.matrix(A)
  n <- nrow(A)

  Q <- matrix(0, n, n)
  R <- matrix(0, n, n)

  for (j in seq_len(n)) {
    v <- A[, j]

    if (j > 1) {
      for (i in seq_len(j - 1)) {
        R[i, j] <- sum(Q[, i] * v)
        v <- v - R[i, j] * Q[, i]
      }
    }

    R[j, j] <- sqrt(sum(v^2))
    if (R[j, j] < tol) stop("Gram-Schmidt QR failed: dependent/near-dependent columns.")

    Q[, j] <- v / R[j, j]
  }

  list(Q = Q, R = R)
}

qr_iteration <- function(A, qr_fun, max_iter = 5000, tol = 1e-12) {
  n <- nrow(A)
  Ak <- A
  Qtotal <- diag(n)

  offdiag_norm <- function(M) {
    M2 <- M
    diag(M2) <- 0
    sqrt(sum(M2^2))
  }

  for (k in 1:max_iter) {
    qr <- qr_fun(Ak)
    Q <- qr$Q; R <- qr$R
    Ak1 <- R %*% Q
    Qtotal <- Qtotal %*% Q

    if (offdiag_norm(Ak1) < tol) {
      Ak <- Ak1
      break
    }
    Ak <- Ak1
  }

  list(values = diag(Ak), vectors = Qtotal)
}

############################################################
# 6) METHOD 5: Householder QR + QR Iteration
############################################################
householder_qr <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  Q <- diag(n)
  R <- A

  for (k in 1:(n - 1)) {
    x <- R[k:n, k]
    e1 <- rep(0, length(x)); e1[1] <- 1

    alpha <- -sign(x[1]) * norm2(x)
    v <- x - alpha * e1
    if (norm2(v) < 1e-15) next
    v <- v / norm2(v)

    Hsmall <- diag(length(x)) - 2 * (v %*% t(v))
    H <- diag(n)
    H[k:n, k:n] <- Hsmall

    R <- H %*% R
    Q <- Q %*% t(H)
  }

  list(Q = Q, R = R)
}

############################################################
# 7) Run all methods and store PC1
############################################################
results <- list()

# Power
pow <- power_iteration(S)
v_pow <- align_sign(pow$vector, v1_true)
results$Power <- list(v1 = v_pow, lam1 = pow$value)

# Deflation
defl <- deflation_power(S)
v_defl <- align_sign(defl$v1, v1_true)
results$Deflation <- list(v1 = v_defl, lam1 = defl$lam1, lam2 = defl$lam2)

# Jacobi
jac <- jacobi_eigen(S)
i1 <- which.max(jac$values)
v_jac <- align_sign(jac$vectors[, i1], v1_true)
results$Jacobi <- list(v1 = v_jac, lam1 = jac$values[i1])

# QR (Gram-Schmidt)
qr_gs_e <- qr_iteration(S, qr_gs)
i1 <- which.max(qr_gs_e$values)
v_qr <- align_sign(qr_gs_e$vectors[, i1], v1_true)
results$QR_GS <- list(v1 = v_qr, lam1 = qr_gs_e$values[i1])

# QR (Householder)
qr_hh_e <- qr_iteration(S, householder_qr)
i1 <- which.max(qr_hh_e$values)
v_hh <- align_sign(qr_hh_e$vectors[, i1], v1_true)
results$QR_Householder <- list(v1 = v_hh, lam1 = qr_hh_e$values[i1])

############################################################
# 8) Display results: loadings + reconstruction error + examples
############################################################
cat("=== Reference (eigen) ===\n")
cat("lambda1_true =", lambda1_true, "\n")
cat("v1_true (PC1 loadings):\n"); print(v1_true)

cat("\n=== Methods (PC1) ===\n")
for (name in names(results)) {
  v1 <- results[[name]]$v1
  lam1 <- results[[name]]$lam1

  rec <- reconstruct_pc1(Xs, v1)
  err <- rmse(Xs, rec$Xhat)

  cos_sim <- abs(sum(v1 * v1_true))  # 1 = same direction
  cat("\n---", name, "---\n")
  cat("lambda1_est =", lam1, "\n")
  cat("v1_est (loadings):\n"); print(v1)
  cat("cosine similarity to true PC1 =", cos_sim, "\n")
  cat("Reconstruction RMSE (PC1 only) =", err, "\n")
}

# Show Real vs Interpolated (first 6 rows) for each method
show_real_vs_hat <- function(method_name) {
  v1 <- results[[method_name]]$v1
  Xhat <- reconstruct_pc1(Xs, v1)$Xhat

  df <- data.frame(
    x1_real = Xs[1:6, 1], x1_hat = Xhat[1:6, 1],
    x2_real = Xs[1:6, 2], x2_hat = Xhat[1:6, 2],
    x3_real = Xs[1:6, 3], x3_hat = Xhat[1:6, 3]
  )

  cat("\nReal vs Interpolated (standardized) -", method_name, "\n")
  print(df)
}

show_real_vs_hat("Power")
show_real_vs_hat("Jacobi")
show_real_vs_hat("QR_Householder")

############################################################
# 9) Plot PCA space (PC1_true vs PC2_true) + reconstructed points
############################################################
scores_true <- Xs %*% cbind(v1_true, v2_true)

plot(scores_true[, 1], scores_true[, 2],
     xlab = "PC1 (true)", ylab = "PC2 (true)",
     main = "PCA space: original points + reconstructed (PC1 only)",
     pch = 16)

# Overlay reconstructed points (each method)
pch_map <- c(Power = 1, Deflation = 2, Jacobi = 3, QR_GS = 4, QR_Householder = 5)

for (name in names(results)) {
  v1 <- results[[name]]$v1
  Xhat <- reconstruct_pc1(Xs, v1)$Xhat
  scores_hat <- Xhat %*% cbind(v1_true, v2_true)  # project to true PCA plane
  points(scores_hat[, 1], scores_hat[, 2], pch = pch_map[name])
}

legend("topright",
       legend = c("Original", names(results)),
       pch = c(16, pch_map[names(results)]),
       bty = "n")
