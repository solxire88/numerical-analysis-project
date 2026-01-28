# Imputation (Interpolation) in R — Direct (Linear) vs Lagrange vs Newton

This mini-project shows a simple way to do **missing value imputation** in a small time series using **interpolation**.

We build a toy “temperature over 7 days” dataset, intentionally remove one value in the middle (a **hole**), then estimate it using:

1. **Direct / naive interpolation** (linear between immediate neighbors)  
2. **Lagrange interpolation polynomial**  
3. **Newton interpolation** (divided differences)

Finally, the script compares **Imputed vs Real** (the value we hid) and plots how each method fills the gap.

---

## Project goal

Given a small set of known points \((x_i, y_i)\), estimate an unknown value at a point \(x^*\):

\[
\text{impute } y(x^*) \text{ from surrounding observations}
\]

---

## Data generation (7-day time series)

- Days are simply:

```r
days <- 1:7
```

- Temperatures follow a smooth pattern + small noise:

```r
temp_true <- 20 + 2*sin(days/2) + rnorm(length(days), 0, 0.15)
```

- We remove the middle value (day 4):

```r
missing_day <- 4
real_hidden <- temp_true[missing_day]  # real value (kept for evaluation)
temp_obs <- temp_true
temp_obs[missing_day] <- NA
```

So `temp_obs` contains a missing value, and `real_hidden` is the correct value we compare against.

---

## Methods implemented

### 1) Direct (Linear) interpolation
Uses the two closest neighbors: **day 3** and **day 5**.

If we have points \((x_0, y_0)\) and \((x_1, y_1)\), then:

\[
y(x) = y_0 + (y_1 - y_0) \frac{x - x_0}{x_1 - x_0}
\]

In code:

```r
impute_direct <- y0 + (y1 - y0) * (x_star - x0) / (x1 - x0)
```

**Pros:** simplest method  
**Cons:** not as accurate if the curve is strongly nonlinear

---

### 2) Lagrange interpolation (degree 3)
We choose 4 surrounding points (to avoid the missing day):

- days: `2, 3, 5, 6`

This fits a **cubic polynomial** (degree 3) through these points.

Lagrange form:

\[
P(x) = \sum_{i=1}^{n} y_i L_i(x)
\]

where:

\[
L_i(x) = \prod_{j \ne i} \frac{x - x_j}{x_i - x_j}
\]

In code, the function `lagrange_eval(x, y, x_star)` evaluates \(P(x^*)\).

---

### 3) Newton interpolation (divided differences)
Newton form:

\[
P(x) = a_0 + a_1(x-x_0) + a_2(x-x_0)(x-x_1) + \dots
\]

The coefficients \(a_k\) come from **divided differences**:

- `newton_divided_diff(x, y)` builds the coefficient vector
- `newton_eval(x, coef, x_star)` evaluates the polynomial at \(x^*\)

**Pros:** efficient to update if you add a new data point  
**Cons:** slightly more code than Lagrange, but very standard in numerical analysis

---

## Output

### Console table
The script prints a table like:

| Method | Imputed | Real | Abs_Error |
|-------|---------|------|----------|
| Direct (Linear) | ... | real_hidden | ... |
| Lagrange (deg 3) | ... | real_hidden | ... |
| Newton (deg 3) | ... | real_hidden | ... |

This shows how close each method is to the true hidden value.

### Plot
The plot displays:

- Observed temperature points (with a missing point at day 4)
- The **real hidden** value at day 4 (marked with a star)
- The **imputed** values (different point symbols)
- Curves showing how each interpolation method fills the gap (between days 2 and 6)

---

## How to run

Save your script as `imputation_interpolation.R` and run:

```bash
Rscript imputation_interpolation.R
```

Or open it in RStudio and click **Run**.

---

## Notes / tips

- The polynomial methods (Lagrange/Newton) can be very accurate **inside** the interval of known points.
- With more points or higher degree, polynomial interpolation can oscillate (Runge phenomenon).  
  Here we keep it stable by using only **4 nearby points**.
- For bigger datasets, people often prefer **splines** (piecewise polynomials), but Lagrange/Newton are perfect for a course demo.


## Author
Made for a school project connecting **interpolation** and **missing data imputation** in R.
