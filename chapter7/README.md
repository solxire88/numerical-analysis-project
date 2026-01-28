# Gradient Descent as an ODE in R — Euler / Taylor(2) / RK4

This mini-project explains a neat viewpoint:

> **Gradient descent is the numerical solution of an ODE.**

Instead of seeing gradient descent only as an optimization loop, we treat it as a **continuous-time dynamical system** and then discretize it using classic **ODE solvers**:

1. **Euler (explicit)** — order 1  
2. **Taylor (order 2)** — uses derivative of the ODE (needs second derivative of the loss)  
3. **Runge–Kutta 4 (RK4)** — order 4  

The script then:
- prints the **final minimum found** (approximately),
- plots the **path taken** by each method on the loss curve,
- and also plots **w(t)** over time.

---

## 1) Loss function (known minimum)

We minimize the simple convex function:

\[
L(w) = (w - 2)^2 + 1
\]

It has a known unique minimum at:

\[
w^* = 2
\]

In the script:

```r
L  <- function(w) (w - 2)^2 + 1
g  <- function(w) 2*(w - 2)     # gradient dL/dw
gp <- function(w) 2             # second derivative of L
```

---

## 2) Gradient descent as an ODE

In continuous time, gradient descent is:

\[
\frac{dw}{dt} = -\frac{dL}{dw} = -g(w)
\]

So the ODE right-hand side is:

```r
f_ode <- function(w) -g(w)
```

---

## 3) Numerical ODE solvers (discrete updates)

We simulate the ODE using step size \(h\). Each method produces a sequence:

\[
w_0, w_1, w_2, \dots
\]

where \(w_k\) approximates \(w(t_k)\) with \(t_k = k h\).

### A) Euler (explicit) — order 1

\[
w_{k+1} = w_k + h f(w_k)
\]

Code:

```r
step_euler <- function(w, h) w + h * f_ode(w)
```

**Pros:** simplest  
**Cons:** less accurate; may require small step size for stability

---

### B) Taylor method (order 2)

For \(w' = f(w)\), the 2nd order Taylor step is:

\[
w_{k+1} = w_k + h f(w_k) + \frac{h^2}{2} f'(w_k) f(w_k)
\]

Here:
- \(f(w) = -g(w)\)
- \(f'(w) = -g'(w)\)
- and \(g'(w)\) is the second derivative of \(L\)

Code:

```r
fp <- -gp(w)
w + h*f + 0.5*h^2*fp*f
```

**Pros:** better accuracy than Euler  
**Cons:** needs second derivative information (or approximation)

---

### C) Runge–Kutta 4 (RK4) — order 4

RK4 uses four slope evaluations:

\[
\begin{aligned}
k_1 &= f(w) \\
k_2 &= f(w + \tfrac{h}{2}k_1) \\
k_3 &= f(w + \tfrac{h}{2}k_2) \\
k_4 &= f(w + hk_3)
\end{aligned}
\]

Update:

\[
w_{k+1} = w + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\]

**Pros:** very accurate for smooth dynamics  
**Cons:** more computations per step

---

## 4) Simulation settings

The script uses:

- Starting point: `w0 = -4`
- Step size: `h = 0.15`
- Steps: `n_steps = 35`

All methods are simulated with the same driver function:

```r
simulate_method(step_fun, w0, h, n_steps)
```

---

## 5) Output

### Console
The script prints:
- the known minimum \(w^*=2\)
- the final estimate from each method:

- `Euler: w_end = ... | L = ...`
- `Taylor2: w_end = ... | L = ...`
- `RK4: w_end = ... | L = ...`

All should be close to \(2\) if the step size is reasonable.

### Plots
1) **Loss curve** \(L(w)\) with the **iterative path** (points) for each method  
2) **w(t)** vs time for each method, plus a horizontal line at \(w^*\)

---

## How to run

Save the script as `gd_ode_methods.R` and run:

```bash
Rscript gd_ode_methods.R
```

Or run inside RStudio.

---

## Notes / troubleshooting

- If the path oscillates or diverges, reduce the step size `h` (e.g., `h = 0.05`).
- Euler is the least accurate; RK4 is the most accurate here.
- Taylor(2) needs the second derivative (here constant `2`, so it’s easy).



## Author
Made for a school project connecting **gradient descent**, **ODEs**, and **numerical integration methods** in R.
