# PCA (Dimensionality Reduction) in R — PC1 via Iterative Eigen Methods

This project demonstrates **PCA (Principal Component Analysis)** on a small synthetic dataset and shows how the **first principal component (PC1)** can be computed using several classic **numerical linear algebra** eigenvalue algorithms.

Instead of calling `prcomp()` to get PCA directly, we:
1. Build the **covariance matrix** \(S\)
2. Compute the dominant eigenvector (PC1) using **iterative methods**
3. Reconstruct (“interpolate”) the data using only PC1 (**rank‑1 approximation**)
4. Compare reconstructions and plot the results in PCA space

---

## 1) Dataset: highly correlated features

We generate **3 correlated variables** \((x_1, x_2, x_3)\) from a latent factor \(z\):

```r
z <- rnorm(n)
X <- cbind(
  x1 = z + 0.20*rnorm(n),
  x2 = 0.9*z + 0.20*rnorm(n),
  x3 = -0.7*z + 0.20*rnorm(n)
)
```

This creates high correlation, so the dataset **compresses well**: most variance is captured by PC1.

### Standardization
We center and scale before PCA:

```r
Xs <- scale(X, center = TRUE, scale = TRUE)
```

### Covariance matrix
PCA is based on the covariance of the standardized data:

```r
S <- cov(Xs)
```

We add a tiny ridge for stability:

```r
S <- S + 1e-10 * diag(nrow(S))
```

---

## 2) PCA as an eigenvalue problem

For covariance matrix \(S\):

- Eigenvalues \(\lambda\) measure variance along components
- Eigenvectors \(v\) are the principal directions

PC1 corresponds to the **largest eigenvalue**:

\[
S v_1 = \lambda_1 v_1, \quad \lambda_1 = \max(\lambda)
\]

The script uses `eigen(S)` only as a **reference (“true” PC1)** for comparison:

```r
eig_ref <- eigen(S)
v1_true <- eig_ref$vectors[, 1]
lambda1_true <- eig_ref$values[1]
```

---

## 3) Methods implemented to compute PC1

The code computes PC1 (dominant eigenvector) with the following methods:

### Method A — Power Iteration (Puissance itérée)
Iteratively applies:

\[
v_{k+1} = \frac{S v_k}{\|S v_k\|}
\]

Converges to the dominant eigenvector when \(\lambda_1\) is separated from others.

**Output:** approximate \(v_1\) and \(\lambda_1\) (Rayleigh quotient).

---

### Method B — Deflation (Déflation)
Uses power iteration to compute \((\lambda_1, v_1)\), then removes that component:

\[
S_1 = S - \lambda_1 \frac{v_1 v_1^T}{v_1^T v_1}
\]

Then power iteration on \(S_1\) estimates the second eigenpair (PC2 demo).

> Deflation is included to show how to compute multiple PCs one-by-one.

---

### Method C — Jacobi rotations (Givens-style) for symmetric matrices
Jacobi repeatedly zeroes the largest off-diagonal element using a rotation:

- Find the largest \(|S_{pq}|\)
- Compute rotation angle \(\theta\)
- Apply a plane rotation to reduce the off-diagonal term

This gradually diagonalizes \(S\), and the accumulated rotation matrix gives eigenvectors.

---

### Method D — QR iteration (Gram–Schmidt QR)
QR iteration does:

\[
S_k = Q_k R_k, \quad S_{k+1} = R_k Q_k
\]

and \(S_k\) converges to an upper triangular matrix whose diagonal approaches eigenvalues.

We build QR using **Gram–Schmidt**, with a safe loop (avoids the R indexing bug for `j=1`):

```r
for (j in seq_len(n)) {
  if (j > 1) { ... }
}
```

---

### Method E — QR iteration (Householder QR)
Same QR iteration idea, but QR is computed using **Householder reflections**, which is typically more numerically stable than classical Gram–Schmidt.

---

## 4) Reconstruction (“Interpolation”) using only PC1

Once a method produces PC1 direction \(v_1\), we project and reconstruct:

### Scores (projection onto PC1)
\[
t = X_s v_1
\]

### Rank‑1 reconstruction
\[
\hat{X} = t v_1^T
\]

In code:

```r
t1 <- Xs %*% v1
Xhat <- t1 %*% t(v1)
```

This is the “compressed” representation using only one dimension.

---

## 5) Displayed results

The script prints:

- Reference \(\lambda_1\), \(v_1\) from `eigen()`
- For each method:
  - Estimated \(\lambda_1\)
  - Estimated \(v_1\)
  - **Cosine similarity** with the reference PC1 (direction match)

\[
\text{cos\_sim} = |v_1^T v_{1,true}|
\]

- Reconstruction error (RMSE):

\[
\text{RMSE} = \sqrt{\frac{1}{np}\sum (X_s - \hat{X})^2}
\]

It also shows “Real vs Interpolated” for the first few rows.

---

## 6) Plot: PCA space + reconstructed points

We plot the original points in the **true PCA plane** (PC1_true vs PC2_true):

```r
scores_true <- Xs %*% cbind(v1_true, v2_true)
```

Then we overlay each method’s reconstructed points projected into the same plane:

```r
scores_hat <- Xhat %*% cbind(v1_true, v2_true)
```

This visualizes how well the rank‑1 approximation preserves structure.

---

## How to run

Save the script as `pca_iterative_methods.R` and run:

```bash
Rscript pca_iterative_methods.R
```

Or run it in RStudio.

---

## Notes / troubleshooting

- **Eigenvector sign is arbitrary**: \(v\) and \(-v\) represent the same component.
  The script aligns signs with the reference (`align_sign()`).
- If Gram–Schmidt QR fails due to near-dependence, the tiny ridge term helps.
- With bigger matrices, prefer Householder QR for stability.



## Author
Made for a school project connecting **PCA** and **iterative eigenvalue algorithms** in R.
