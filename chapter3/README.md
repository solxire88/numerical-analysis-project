# 1D “Gradient Descent” via Root-Finding (Bisection / Fixed Point / Newton) — R

This mini-project demonstrates a common optimization idea in 1D:

> **To minimize a loss function \(L(w)\), we can solve its gradient equation**  
> \(\frac{dL}{dw} = 0\).

Instead of using a full gradient descent loop, we **directly find the stationary point** (where the derivative is zero) using three classic iterative methods:

1. **Dichotomy / Bisection** (bracketing method)
2. **Point fixe / Fixed-point iteration**
3. **Newton–Raphson** (uses derivative information)

Finally, the script plots \(L(w)\) and marks the solutions found by each method.

---

## Model (loss function)

We minimize the function:

\[
L(w) = \tfrac{1}{2} w^2 + \sin(w)
\]

### Gradient (first derivative)
\[
f(w) = \frac{dL}{dw} = w + \cos(w)
\]

We want to find:

\[
f(w) = 0
\]

### Second derivative (derivative of the gradient)
\[
f'(w) = 1 - \sin(w)
\]

This is used by **Newton–Raphson**.

---

## Methods implemented

All three methods aim to solve the same equation:

\[
f(w) = 0
\]

They return:
- `root` : the estimated solution \(w^*\)
- `iters` : number of iterations needed

### 1) Dichotomy (Bisection)

**When to use:**  
When you can find an interval \([a,b]\) such that:

\[
f(a) \cdot f(b) < 0
\]

That means a root lies inside the interval (sign change).

**Update rule:**  
Take the midpoint:

\[
m = \frac{a+b}{2}
\]

and keep the half-interval where the sign change remains.

**Pros:** guaranteed convergence (if bracket exists)  
**Cons:** slower than Newton

In the script:

```r
sol_bis <- dichotomy(f, a = -2, b = 1)
```

---

### 2) Point fixe (Fixed-point iteration)

We rewrite the root equation into:

\[
w = g(w)
\]

A common choice is:

\[
g(w) = w - \alpha f(w)
\]

This looks like a “gradient descent” step with step size \(\alpha\).

**Update rule:**
\[
w_{k+1} = g(w_k)
\]

In the script:

```r
alpha <- 0.005
g <- function(w) w - alpha * f(w)
sol_fp <- point_fixe(g, w0 = 0)
```

**Pros:** simple (no derivative needed)  
**Cons:** converges only if \(g\) is a contraction (step size must be small enough)

---

### 3) Newton–Raphson

Uses both \(f(w)\) and \(f'(w)\):

\[
w_{k+1} = w_k - \frac{f(w_k)}{f'(w_k)}
\]

In the script:

```r
sol_new <- newton_raphson(f, fp, w0 = 0)
```

**Pros:** very fast (quadratic convergence near the solution)  
**Cons:** can fail if \(f'(w)\) is close to 0 or if the initial guess is bad

---

## How to run

Save your script as, for example, `root_minimize.R`, then run:

```bash
Rscript root_minimize.R
```

Or open in RStudio and click **Run**.

---

## Output

### Console output
The script prints the minimizer found by each method and the final loss value:

- `Bisection      w* = ... | iters = ... | L(w*) = ...`
- `Point fixe     w* = ... | iters = ... | L(w*) = ...`
- `Newton-Raphson w* = ... | iters = ... | L(w*) = ...`

Because all methods solve the same equation \(f(w)=0\), the solutions should be very close.

### Plot
The plot shows:
- the curve \(L(w)\)
- three markers for the solutions found by:
  - Bisection
  - Fixed point
  - Newton–Raphson

---

## Notes / troubleshooting

- **Bisection needs a valid bracket**: if you choose \([a,b]\) with no sign change, it stops.
- **Fixed point depends on `alpha`**:
  - if it’s too large, it may diverge
  - smaller values are safer but slower
- **Newton needs `fp(w)` not too small**:
  - the code stops if `abs(fp(w)) < 1e-14`



## Author
Made for a school project linking **optimization** and **iterative numerical methods** in R.
