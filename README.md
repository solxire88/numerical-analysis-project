# Numerical Methods for Machine Learning (R) — Mini Project Collection

This repository is a collection of small **school-friendly** R scripts that connect **numerical analysis** (linear algebra, iterative solvers, interpolation, numerical integration, ODE solvers) with **machine learning concepts** (regression, PCA, ROC/AUC, imputation, optimization).

Each script is **self-contained**:
- generates its own synthetic data,
- implements the numerical method(s) from scratch,
- prints results,
- and produces at least one plot.

---

## Table of Contents

1. [Requirements](#requirements)  
3. [Projects Included](#projects-included)  
   - [A) Polynomial Regression — Direct Solvers (Gauss / LU / Cholesky)](#a-polynomial-regression--direct-solvers-gauss--lu--cholesky)  
   - [B) Polynomial Regression — Iterative Solvers (Jacobi / Gauss–Seidel / SOR)](#b-polynomial-regression--iterative-solvers-jacobi--gaussseidel--sor)  
   - [C) 1D Minimization via Root-Finding (Bisection / Fixed-Point / Newton)](#c-1d-minimization-via-root-finding-bisection--fixed-point--newton)  
   - [D) PCA — PC1 via Iterative Eigen Methods (Power / Deflation / Jacobi / QR / Householder)](#d-pca--pc1-via-iterative-eigen-methods-power--deflation--jacobi--qr--householder)  
   - [E) ROC / AUC — Numerical Integration (Midpoint / Trapezoid / Simpson)](#e-roc--auc--numerical-integration-midpoint--trapezoid--simpson)  
   - [F) Time-Series Imputation — Interpolation (Direct / Lagrange / Newton)](#f-time-series-imputation--interpolation-direct--lagrange--newton)  
   - [G) Gradient Descent as an ODE (Euler / Taylor(2) / RK4)](#g-gradient-descent-as-an-ode-euler--taylor2--rk4)   

---

## Requirements

- **R** (any recent version should work)
- No external libraries required (base R only)

Optional (recommended):
- **RStudio** to run scripts interactively



## Projects Included



### A) Polynomial Regression — Direct Solvers (Gauss / LU / Cholesky)

**Goal:** Train a polynomial regression model with 4 coefficients by solving the normal equation:

\[
(X^T X) w = X^T y
\]

**Numerical methods implemented:**
- Gaussian elimination (with partial pivoting)
- LU decomposition (Doolittle + partial pivoting)
- Cholesky decomposition \(A = L L^T\) (with a tiny ridge term for stability)

**Output:**
- printed coefficient table (true vs estimated)
- plot of data + fitted curve(s)

---

### B) Polynomial Regression — Iterative Solvers (Jacobi / Gauss–Seidel / SOR)

**Goal:** Solve the same normal equation system \(A w = b\), but using iterative methods.

**Numerical methods implemented:**
- Jacobi method
- Gauss–Seidel method
- SOR (Successive Over-Relaxation)

**Stability note:** a tiny ridge term is added to make \(A\) better conditioned and improve convergence.


**Output:**
- coefficient table
- iteration counts + residual norm
- plot of data + fitted curves

---

### C) 1D Minimization via Root-Finding (Bisection / Fixed-Point / Newton)

**Goal:** Minimize a 1D loss \(L(w)\) by solving the gradient equation:

\[
\frac{dL}{dw} = 0
\]

**Model used:**
\[
L(w) = 0.5 w^2 + \sin(w)
\]

**Numerical methods implemented:**
- Dichotomy / Bisection on \(f(w)=0\)
- Fixed-point iteration \(w = g(w)\)
- Newton–Raphson


**Output:**
- solutions \(w^*\) and iteration counts
- plot of \(L(w)\) with markers at the found minima

---

### D) PCA — PC1 via Iterative Eigen Methods (Power / Deflation / Jacobi / QR / Householder)

**Goal:** Perform PCA by computing the dominant eigenvector of the covariance matrix \(S\).

**Pipeline:**
1. Generate correlated features
2. Standardize
3. Build covariance matrix \(S = cov(X)\)
4. Compute **PC1** with multiple eigen algorithms
5. Reconstruct data using only PC1 (rank-1 approximation)

**Methods implemented:**
- Power iteration
- Deflation (PC2 demo)
- Jacobi rotations
- QR iteration with Gram–Schmidt QR
- QR iteration with Householder QR


**Output:**
- eigenvalue/eigenvector estimates vs reference
- reconstruction RMSE
- PCA-space plot: original points + reconstructed points

---

### E) ROC / AUC — Numerical Integration (Midpoint / Trapezoid / Simpson)

**Goal:** Compute AUC (area under ROC curve) by numerical integration.

**ROC curve used:**
\[
TPR = FPR^{\alpha},\quad 0<\alpha<1
\]
This has an exact AUC:
\[
AUC_{true} = \int_0^1 x^{\alpha} dx = \frac{1}{\alpha+1}
\]

**Methods implemented:**
- Composite midpoint (rectangles) quadrature
- Composite trapezoidal rule
- Composite Simpson rule


**Output:**
- table: exact AUC vs approximate AUC and absolute errors
- ROC plot with shaded area

---

### F) Time-Series Imputation — Interpolation (Direct / Lagrange / Newton)

**Goal:** Fill a missing value in a small 7-day temperature series.

**Imputation methods:**
- Direct (linear) interpolation between neighbors
- Lagrange polynomial interpolation (degree 3)
- Newton interpolation (divided differences, degree 3)


**Output:**
- table: imputed vs real hidden value + absolute error
- plot: series + interpolation curves and imputed point

---

### G) Gradient Descent as an ODE (Euler / Taylor(2) / RK4)

**Goal:** Show gradient descent as a continuous-time ODE:

\[
\frac{dw}{dt} = -\nabla L(w)
\]

Then discretize it using ODE methods.

**Loss function (known minimum):**
\[
L(w) = (w-2)^2 + 1\quad \Rightarrow\quad w^*=2
\]

**ODE solvers implemented:**
- Euler (order 1)
- Taylor order 2 (needs second derivative)
- Runge–Kutta 4 (RK4)


**Output:**
- final estimated minimum for each method
- plot: loss curve with the iteration paths
- plot: w(t) over time


