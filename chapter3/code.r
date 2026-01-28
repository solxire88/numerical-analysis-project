# ============================================================
# Simple "Gradient Descent" idea in 1D (R)
# We minimize a loss L(w) by solving its gradient = 0:
#   f(w) = dL/dw = 0
#
# We find the solution w* using 3 iterative methods:
#   1) Dichotomy (Bisection) on f(w)=0
#   2) Point fixe (Fixed-point iteration) on w = g(w)
#   3) Newton–Raphson on f(w)=0
#
# Then we plot L(w) and mark the 3 solutions.
# ============================================================

# ----------------------------
# 1) Define a simple model (loss) and its derivatives
#    New Model: L(w) = 0.5*w^2 + sin(w)
# ----------------------------
L  <- function(w) 0.5 * w^2 + sin(w)

# Gradient: f(w) = w + cos(w)
# We want to find w where f(w) = 0
f  <- function(w) w + cos(w)

# Derivative of gradient: f'(w) = 1 - sin(w)
# (Needed for Newton-Raphson)
fp <- function(w) 1 - sin(w)


# ----------------------------
# 2) Method 1: Dichotomy (Bisection) for f(w)=0
#    Needs an interval [a,b] where f(a)*f(b) < 0
# ----------------------------
dichotomy <- function(f, a, b, tol = 1e-10, max_iter = 200) {
  fa <- f(a); fb <- f(b)
  if (fa * fb > 0) stop("Bisection needs f(a)*f(b) < 0 (root must be bracketed).")

  for (k in 1:max_iter) {
    m <- (a + b) / 2
    fm <- f(m)

    if (abs(fm) < tol || (b - a)/2 < tol) {
      return(list(root = m, iters = k))
    }

    if (fa * fm < 0) {
      b <- m; fb <- fm
    } else {
      a <- m; fa <- fm
    }
  }
  list(root = (a + b)/2, iters = max_iter)
}


# ----------------------------
# 3) Method 2: Point fixe (Fixed-point)
#    Build g(w) so that f(w)=0 <=> g(w)=w
#    A common choice: g(w) = w - alpha*f(w)
#    If alpha is small enough, it converges.
# ----------------------------
point_fixe <- function(g, w0, tol = 1e-10, max_iter = 5000) {
  w <- w0
  for (k in 1:max_iter) {
    w_next <- g(w)
    if (abs(w_next - w) < tol) {
      return(list(root = w_next, iters = k))
    }
    w <- w_next
  }
  list(root = w, iters = max_iter)
}

alpha <- 0.005
g <- function(w) w - alpha * f(w)


# ----------------------------
# 4) Method 3: Newton–Raphson for f(w)=0
#    w_{n+1} = w_n - f(w_n)/f'(w_n)
# ----------------------------
newton_raphson <- function(f, fp, w0, tol = 1e-10, max_iter = 100) {
  w <- w0
  for (k in 1:max_iter) {
    if (abs(fp(w)) < 1e-14) stop("Newton stopped: derivative too small.")
    w_next <- w - f(w)/fp(w)
    if (abs(w_next - w) < tol) {
      return(list(root = w_next, iters = k))
    }
    w <- w_next
  }
  list(root = w, iters = max_iter)
}


# ============================================================
# 5) Run the 3 methods on the gradient equation f(w)=0
# ============================================================

# Bisection bracket (here it works because f(0)<0 and f(4)>0)
sol_bis <- dichotomy(f, a = -2, b = 1)

# Fixed point needs an initial guess
sol_fp  <- point_fixe(g, w0 = 0)

# Newton-Raphson initial guess
sol_new <- newton_raphson(f, fp, w0 = 0)

w_bis <- sol_bis$root
w_fix <- sol_fp$root
w_new <- sol_new$root

cat("Solutions (minimizer w* where dL/dw = 0)\n")
cat("Bisection      w* =", w_bis, " | iters =", sol_bis$iters, " | L(w*) =", L(w_bis), "\n")
cat("Point fixe     w* =", w_fix, " | iters =", sol_fp$iters,  " | L(w*) =", L(w_fix), "\n")
cat("Newton-Raphson w* =", w_new, " | iters =", sol_new$iters, " | L(w*) =", L(w_new), "\n")


# ============================================================
# 6) Plot the loss function and mark the 3 solutions
# ============================================================

wg <- seq(-1, 5, length.out = 500)
Lg <- L(wg)

plot(wg, Lg, type = "l",
     main = "Minimize L(w) by solving dL/dw = 0 (3 methods)",
     xlab = "w", ylab = "L(w)")

# Mark each computed solution on the curve
points(w_bis, L(w_bis), pch = 19)
points(w_fix, L(w_fix), pch = 17)
points(w_new, L(w_new), pch = 15)

legend("topright",
       legend = c("Bisection", "Point fixe", "Newton-Raphson"),
       pch = c(19, 17, 15),
       bty = "n")
