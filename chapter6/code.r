############################################################
# AUC (Area Under ROC Curve) with 3 numerical integration methods
# - Create sample ROC data (FPR vs TPR) from a known curve
# - Compute AUC using:
#     1) Formule de quadrature (ici: rectangles / midpoint composite)
#     2) Trapezoidal (composite)
#     3) Simpson (composite)
# - Compare with the "real" AUC (analytical integral)
# - Plot ROC curve + shade the AUC area
############################################################

set.seed(42)

#############################
# 1) Create a sample ROC curve (highly "good" classifier shape)
#    We choose: TPR = (FPR)^alpha, with 0 < alpha < 1 => concave ROC.
#    Real AUC = ∫0^1 x^alpha dx = 1/(alpha+1)
#############################
alpha <- 0.35
fpr <- seq(0, 1, length.out = 101)     # 100 intervals (good for Simpson)
tpr <- fpr^alpha                       # ROC points

auc_true <- 1 / (alpha + 1)            # "Real" AUC (exact)

#############################
# 2) AUC with numerical methods (on discrete points)
#############################

# --- (1) Quadrature formula (Rectangles / Midpoint composite)
# We approximate TPR at each interval midpoint using the average of endpoints
auc_midpoint <- function(x, y) {
  n <- length(x)
  h <- x[2] - x[1]                     # assume equally spaced
  y_mid <- (y[1:(n-1)] + y[2:n]) / 2   # midpoint value by linear approximation
  sum(h * y_mid)
}

# --- (2) Trapezoidal composite
auc_trapezoid <- function(x, y) {
  n <- length(x)
  h <- x[2] - x[1]
  h * ( (y[1] + y[n]) / 2 + sum(y[2:(n-1)]) )
}

# --- (3) Simpson composite
# Needs an even number of intervals => (n-1) must be even
auc_simpson <- function(x, y) {
  n <- length(x)
  m <- n - 1                           # number of intervals
  if (m %% 2 != 0) stop("Simpson: need an even number of intervals (n-1 must be even).")
  h <- x[2] - x[1]

  # indices: 1..n
  # odd interior points: 2,4,6,...,n-1  (in R indexing these are even indices)
  # even interior points: 3,5,7,...,n-2 (in R indexing these are odd indices)
  odd_idx  <- seq(2, n-1, by = 2)
  even_idx <- seq(3, n-2, by = 2)

  (h/3) * ( y[1] + y[n] + 4*sum(y[odd_idx]) + 2*sum(y[even_idx]) )
}

# Compute AUCs
auc_q  <- auc_midpoint(fpr, tpr)
auc_t  <- auc_trapezoid(fpr, tpr)
auc_s  <- auc_simpson(fpr, tpr)

#############################
# 3) Display results
#############################
results <- data.frame(
  method = c("Real (exact)", "Quadrature (midpoint)", "Trapezoidal", "Simpson"),
  auc    = c(auc_true, auc_q, auc_t, auc_s),
  abs_error = c(0, abs(auc_q - auc_true), abs(auc_t - auc_true), abs(auc_s - auc_true))
)

print(results)

#############################
# 4) Plot ROC curve + shade AUC
#############################
plot(fpr, tpr, type = "l", lwd = 2,
     xlab = "False Positive Rate (FPR)",
     ylab = "True Positive Rate (TPR)",
     main = "ROC Curve with Shaded AUC")

# Shade area under ROC curve
polygon(
  c(fpr, rev(fpr)),
  c(tpr, rep(0, length(tpr))),
  col = rgb(0.2, 0.6, 0.9, 0.25),      # semi-transparent fill
  border = NA
)

# Draw curve again on top
lines(fpr, tpr, lwd = 2)

# Add diagonal baseline (random classifier)
abline(0, 1, lty = 2)

legend("bottomright",
       legend = c(
         paste0("AUC true = ", round(auc_true, 4)),
         paste0("Midpoint = ", round(auc_q, 4)),
         paste0("Trapezoid = ", round(auc_t, 4)),
         paste0("Simpson = ", round(auc_s, 4))
       ),
       bty = "n")
