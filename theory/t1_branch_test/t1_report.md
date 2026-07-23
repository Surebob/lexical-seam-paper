# T1 Branch-Structure Test

The shifted multilingual winner passes the three requested Bregman boundary checks. For `f(x) = x^x - sqrt(x) - 0.5(x - 1)`, symbolic evaluation gives `f(1) = 0` and `f'(1) = 0`, matching the handwritten derivative `x^x(1 + log x) - 1/(2 sqrt x) - 0.5`. The second derivative simplifies to `x^x*((log(x)+1)^2 + 1/x) + 1/(4*x^(3/2))`, which is strictly positive for every `x > 0`, hence on `[0.05, 1.0]`.

The 1000-point numerical grid agrees with the symbolic result: `|f(1)| = 0`, `|f'(1)| = 0`, and every sampled `f''(x)` is positive, with minimum `1.86124` at `x = 0.713764`. As a bonus diagnostic, the shifted branch does not reduce cleanly to the IS, exponential, or Euclidean generators under a simple centered-shape comparison; it is best treated as its own convex generator branch rather than one of those known forms.
