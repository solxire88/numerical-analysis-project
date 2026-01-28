# Polynomial Linear Regression (4 coefficients) trained with Jacobi / Gauss–Seidel / SOR (R)

This mini-project demonstrates how **iterative methods for linear systems** can be used to train a simple **machine learning regression model** *without* calling `lm()` or directly inverting matrices.

We fit a polynomial regression model with **4 coefficients**:

\[
y = w_0 + w_1 x + w_2 x^2 + w_3 x^3 + \text{noise}
\]

Training is done by solving the **Normal Equation**:

\[
(X^T X) w = X^T y
\]

Let:

- \(A = X^T X\)
- \(b = X^T y\)

Then training becomes solving the linear system:

\[
A w = b
\]

We solve it using three classic iterative methods:

1. **Jacobi**
2. **Gauss–Seidel**
3. **SOR (Successive Over-Relaxation)**

---

## What the script does

### 1) Synthetic dataset (simple + controlled)
- Generates `n = 140` random values `x ∈ [-2, 2]`
- Uses a known “true” coefficient vector:

```r
w_true <- c(2.0, -1.5, 0.7, 0.2)
```

- Builds the design matrix:

```r
X <- cbind(1, x, x^2, x^3)
```

- Generates noisy targets:

```r
y <- as.numeric(X %*% w_true + noise)
```

This gives you a dataset where you already know what the coefficients **should** be.

---

### 2) Normal equation system
The code forms:

```r
A <- t(X) %*% X
b <- t(X) %*% y
```

So we solve:

\[
A w = b
\]

#### Why we add a small ridge term
Iterative methods converge more reliably when the matrix is well-conditioned.
Since \(A = X^T X\) is symmetric and (often) positive semi-definite, we add a tiny ridge:

```r
ridge <- 1e-8
A <- A + ridge * diag(nrow(A))
```

This makes `A` **strictly positive definite** in practice and helps convergence.

---

## Iterative methods used

All methods return:
- `w` : estimated coefficients
- `iters` : number of iterations
- `residual` : \(\|Aw - b\|\) at the end (how close we are)

### 1) Jacobi
Update rule (component-wise):

\[
w_i^{(k+1)} = \frac{b_i - \sum_{j\ne i} A_{ij} w_j^{(k)}}{A_{ii}}
\]

- Uses only the previous iterate values (slower, but simple).

### 2) Gauss–Seidel
Update rule:

\[
w_i^{(k+1)} = \frac{b_i - \sum_{j<i} A_{ij} w_j^{(k+1)} - \sum_{j>i} A_{ij} w_j^{(k)}}{A_{ii}}
\]

- Uses newly updated values immediately (usually converges faster than Jacobi).

### 3) SOR (Relaxation)
SOR is Gauss–Seidel + a relaxation parameter \(\omega\):

\[
w_i \leftarrow (1-\omega) w_i^{old} + \omega\, w_i^{GS}
\]

- `omega = 1` → Gauss–Seidel
- `1 < omega < 2` → often faster (over-relaxation)
- `0 < omega < 1` → more stable (under-relaxation)

In the script:

```r
omega = 1.2
```

---

## How to run

If your file is named `main.R`:

```bash
Rscript main.R
```

Or open in RStudio and click **Run**.

---

## Output

### Console output
1) A coefficient table comparing:
- true coefficients vs estimated coefficients (Jacobi / GS / SOR)

2) Iterations and residuals, for example:

- `Jacobi : ... iters | residual = ...`
- `Gauss-Seidel : ...`
- `SOR : ...`

### Plot
A scatter plot of the data plus three fitted curves:

- Jacobi
- Gauss–Seidel
- SOR

They should almost overlap if the methods converge well.

---

## Notes / troubleshooting

- If Jacobi is slow, try:
  - increasing `max_iter`
  - slightly increasing `tol` (less strict)
- If SOR diverges, reduce `omega` closer to 1 (e.g., `omega = 1.05`).
- The ridge term `1e-8` is intentionally tiny; it improves conditioning without changing results much.


## Author
Made for a school project linking **iterative numerical methods** and **machine learning**.
