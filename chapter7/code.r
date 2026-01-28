############################################################
# Gradient Descent as an ODE + 3 iterative ODE solvers in R
#
# Idea:
#   Gradient descent in continuous time is the ODE:
#       dw/dt = - dL/dw
#
# If we discretize time with step h, we get iterative methods.
# We will solve the ODE with 3 methods:
#   1) Euler (explicit)
#   2) Taylor order 2 (uses second derivative)
#   3) Runge–Kutta order 4 (RK4)
#
# Then:
# - Display the final minimum found (w_end and L(w_end))
# - Plot the path w(t) taken by each method over the loss curve
############################################################

set.seed(42)

#############################
# 1) Choose a simple loss curve with known minimum
#
# L(w) = (w - 2)^2 + 1
# Minimum at w* = 2
#############################
L  <- function(w) (w - 2)^2 + 1          # loss
g  <- function(w) 2*(w - 2)              # gradient dL/dw
gp <- function(w) 2                      # derivative of gradient (second derivative of L)

w_star <- 2

#############################
# 2) Gradient descent ODE
#   dw/dt = -g(w)
#############################
f_ode <- function(w) -g(w)   # ODE RHS

############################################################
# 3) ODE numerical methods (one step)
############################################################

# 3.1 Euler (order 1)
# w_{k+1} = w_k + h * f(w_k)
step_euler <- function(w, h) {
  w + h * f_ode(w)
}

# 3.2 Taylor order 2
# For ODE w' = f(w), the 2nd order Taylor step is:
# w_{k+1} = w + h f(w) + (h^2/2) f'(w) f(w)
#
# Here f(w) = -g(w)
# f'(w) = d/dw (-g(w)) = -g'(w)
step_taylor2 <- function(w, h) {
  f  <- f_ode(w)
  fp <- -gp(w)             # derivative of f(w)
  w + h * f + 0.5 * h^2 * fp * f
}

# 3.3 Runge–Kutta 4 (order 4)
# Standard RK4:
# k1 = f(w)
# k2 = f(w + h/2 k1)
# k3 = f(w + h/2 k2)
# k4 = f(w + h k3)
# w_{k+1} = w + (h/6)(k1 + 2k2 + 2k3 + k4)
step_rk4 <- function(w, h) {
  k1 <- f_ode(w)
  k2 <- f_ode(w + 0.5*h*k1)
  k3 <- f_ode(w + 0.5*h*k2)
  k4 <- f_ode(w + h*k3)
  w + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
}

############################################################
# 4) Simulate the path for each method
############################################################
simulate_method <- function(step_fun, w0, h, n_steps) {
  w_path <- numeric(n_steps + 1)
  w_path[1] <- w0
  for (k in 1:n_steps) {
    w_path[k + 1] <- step_fun(w_path[k], h)
  }
  w_path
}

# Starting point (far from minimum)
w0 <- -4

# Time step and number of steps
h <- 0.15
n_steps <- 35
t_grid <- seq(0, n_steps*h, by = h)

# Paths
w_euler  <- simulate_method(step_euler,  w0, h, n_steps)
w_taylor <- simulate_method(step_taylor2, w0, h, n_steps)
w_rk4    <- simulate_method(step_rk4,    w0, h, n_steps)

############################################################
# 5) Display final results
############################################################
cat("Known minimum: w* =", w_star, " | L(w*) =", L(w_star), "\n\n")

cat("Final (after", n_steps, "steps, h =", h, ")\n")
cat("Euler:    w_end =", round(tail(w_euler, 1), 6),  " | L =", round(L(tail(w_euler, 1)), 6), "\n")
cat("Taylor2:  w_end =", round(tail(w_taylor,1), 6),  " | L =", round(L(tail(w_taylor,1)), 6), "\n")
cat("RK4:      w_end =", round(tail(w_rk4, 1), 6),    " | L =", round(L(tail(w_rk4, 1)), 6), "\n")

############################################################
# 6) Plot: loss curve + the path taken by each method
############################################################

# Loss curve for plotting
w_grid <- seq(min(c(w_euler, w_taylor, w_rk4)) - 0.5,
              max(c(w_euler, w_taylor, w_rk4)) + 0.5,
              length.out = 400)

plot(w_grid, L(w_grid), type = "l", lwd = 2,
     main = "Gradient Descent Path via ODE Solvers (Euler / Taylor2 / RK4)",
     xlab = "w", ylab = "L(w)")

# Mark true minimum
points(w_star, L(w_star), pch = 19)
text(w_star, L(w_star), labels = "true min", pos = 3)

# Plot paths on the loss curve (each point is an iteration)
points(w_euler,  L(w_euler),  pch = 1)
points(w_taylor, L(w_taylor), pch = 2)
points(w_rk4,    L(w_rk4),    pch = 3)

legend("topright",
       legend = c("Euler", "Taylor order 2", "RK4", "True minimum"),
       pch = c(1, 2, 3, 19),
       bty = "n")

# Optional: also plot w(t) vs time (path over time)
plot(t_grid, w_euler, type = "l", lwd = 2,
     main = "w(t) path for each method",
     xlab = "t", ylab = "w")
lines(t_grid, w_taylor, lwd = 2, lty = 2)
lines(t_grid, w_rk4, lwd = 2, lty = 3)
abline(h = w_star, lty = 2)
legend("topright",
       legend = c("Euler", "Taylor2", "RK4", "w*"),
       lty = c(1, 2, 3, 2),
       lwd = c(2, 2, 2, 1),
       bty = "n")
