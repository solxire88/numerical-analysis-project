# AUC (Area Under the ROC Curve) in R — Midpoint / Trapezoid / Simpson

This project demonstrates how to compute **AUC (Area Under the ROC Curve)** using classic **numerical integration** methods from calculus/numerical analysis.

We:
1. Generate a synthetic ROC curve \((\text{FPR}, \text{TPR})\) from a known function
2. Compute AUC using **three** numerical methods:
   - Composite **midpoint (rectangles)** quadrature
   - Composite **trapezoidal** rule
   - Composite **Simpson** rule
3. Compare numerical results to the **true AUC** (analytical integral)
4. Plot the ROC curve and **shade the AUC area**

---

## What is AUC?

An ROC curve plots:
- **FPR** (False Positive Rate) on the x-axis
- **TPR** (True Positive Rate) on the y-axis

AUC is the area under this curve:

\[
\text{AUC} = \int_0^1 \text{TPR}(\text{FPR})\, d(\text{FPR})
\]

- AUC = 0.5 means random classifier
- AUC closer to 1 means better classifier

---

## 1) Synthetic ROC curve with known “real” AUC

We create ROC points from a concave curve:

\[
\text{TPR} = (\text{FPR})^{\alpha}, \quad 0 < \alpha < 1
\]

This produces a “good classifier” shape (above the diagonal).

### True AUC (analytical integral)

\[
\text{AUC}_{true} = \int_0^1 x^{\alpha} dx = \frac{1}{\alpha + 1}
\]

In the script:

```r
alpha <- 0.35
fpr <- seq(0, 1, length.out = 101)
tpr <- fpr^alpha
auc_true <- 1 / (alpha + 1)
```

We use 101 points (100 intervals), which is convenient for Simpson’s rule (needs an even number of intervals).

---

## 2) Numerical integration methods

All methods approximate:

\[
\int_0^1 \text{TPR}(x)\,dx
\]

based on discrete samples \((x_i, y_i)\), where:
- \(x_i = \text{FPR}_i\)
- \(y_i = \text{TPR}_i\)

### A) Quadrature (Rectangles / Midpoint composite)

We use a composite midpoint rule:
- Step size: \(h\)
- Approximate midpoint values using the average of endpoints:

\[
y_{mid,i} \approx \frac{y_i + y_{i+1}}{2}
\]

Then:

\[
\text{AUC} \approx \sum_i h\, y_{mid,i}
\]

In code:

```r
y_mid <- (y[1:(n-1)] + y[2:n]) / 2
sum(h * y_mid)
```

---

### B) Trapezoidal rule (composite)

For equally spaced samples:

\[
\text{AUC} \approx h\left(\frac{y_0 + y_n}{2} + \sum_{i=1}^{n-1} y_i\right)
\]

In code:

```r
h * ( (y[1] + y[n]) / 2 + sum(y[2:(n-1)]) )
```

---

### C) Simpson rule (composite)

Simpson’s rule requires an **even number of intervals** (here 100).

Formula:

\[
\text{AUC} \approx \frac{h}{3}
\left(
y_0 + y_n
+ 4\sum y_{odd}
+ 2\sum y_{even}
\right)
\]

In the script, indexes are chosen carefully for R’s 1-based indexing:

```r
odd_idx  <- seq(2, n-1, by = 2)
even_idx <- seq(3, n-2, by = 2)
```

---

## 3) Output

The script prints a table:

- **Real (exact)** AUC
- AUC by:
  - Midpoint quadrature
  - Trapezoidal
  - Simpson
- Absolute error for each numerical method

Example columns:

- `method`
- `auc`
- `abs_error`

---

## 4) Plot

The plot shows:
- The ROC curve
- A shaded polygon under the curve (AUC area)
- The diagonal baseline \(y=x\) (random classifier)
- A legend with the AUC values for each method

---

## How to run

Save the script as `auc_integration.R` and run:

```bash
Rscript auc_integration.R
```

Or run it inside RStudio.

---

## Notes / troubleshooting

- **Simpson requires an even number of intervals**:
  - If you change the number of points, ensure `(length(fpr) - 1)` is even.
- The midpoint method here approximates midpoints using endpoint averages; this is simple and works well on smooth curves.
- On real ROC curves from real classifier scores, points might not be evenly spaced. In that case:
  - Trapezoidal rule is the most common practical approach.



## Author
Made for a school project linking **ROC/AUC** and **numerical integration** methods in R.
