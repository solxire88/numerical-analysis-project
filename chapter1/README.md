# Polynomial Linear Regression (4 coefficients) solved with Gauss / LU / Cholesky (R)

This mini-project shows **how classic numerical linear algebra methods** can be used to **train a machine learning model**.

We fit a **polynomial linear regression** model with **4 coefficients**:

\[
y = w_0 + w_1 x + w_2 x^2 + w_3 x^3 + \text{noise}
\]

Instead of using `lm()` or matrix inversion, we compute the weights by solving the **Normal Equation**:

\[
(X^T X)\, w = X^T y
\]

where the design matrix is:

\[
X = [1, x, x^2, x^3]
\]

We solve the system using three different methods:

1. **Gaussian Elimination (Gauss)** with partial pivoting  
2. **LU Decomposition** (Doolittle) with partial pivoting  
3. **Cholesky Decomposition** \(A = L L^T\) (with a tiny ridge for stability)

---

## What the script does

### 1) Generate synthetic training data
- Generates `n = 120` random x-values in `[-2, 2]`
- Defines true coefficients:

```r
w_true <- c(2.0, -1.5, 0.7, 0.2)
```

- Builds the feature matrix:

```r
X <- cbind(1, x, x^2, x^3)
```

- Produces targets with noise:

```r
y <- as.numeric(X %*% w_true + noise)
```

### 2) Build the Normal Equation system
Creates:

- \(A = X^T X\)
- \(b = X^T y\)

So the problem becomes:

\[
A w = b
\]

### 3) Solve \(A w = b\) using 3 methods

#### (a) Gauss (Gaussian Elimination)
- Uses an **augmented matrix** `[A | b]`
- Performs:
  - **partial pivoting** (swap rows to avoid numerical issues)
  - **forward elimination**
  - **back substitution**

#### (b) LU Decomposition (Doolittle)
- Factorizes:

\[
P A = L U
\]

- Solves in two steps:
  1. **Forward substitution**: \(L y = P b\)
  2. **Back substitution**: \(U w = y\)

#### (c) Cholesky
For symmetric positive definite matrices:

\[
A = L L^T
\]

Here `A = X^T X` is symmetric. To ensure it’s strictly positive definite (especially with noisy / correlated data),
we add a tiny **ridge term**:

```r
A_ridge <- A + 1e-10 * diag(nrow(A))
```

Then we solve:
1. \(L z = b\)
2. \(L^T w = z\)

---

## How to run

In a terminal:

```bash
Rscript main.R
```

Or inside R / RStudio, just run the script.

---

## Output

### Console output
- A table showing the **true coefficients** vs estimated coefficients from:
  - Gauss
  - LU
  - Cholesky

### Plot
- Scatter of generated points `(x, y)`
- Three fitted polynomial curves:
  - Gauss
  - LU
  - Cholesky

They should nearly overlap.

---

## Notes / Tips

- **All three methods should give almost the same weights**, because they solve the same linear system.
- Cholesky is typically fastest **when the matrix is SPD**.
- If you set the ridge term to `0`, Cholesky may fail if `A` is not positive definite due to numerical issues.

