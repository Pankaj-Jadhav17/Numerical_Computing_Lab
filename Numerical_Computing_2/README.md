**Project Overview**

This repository implements numerical-computing tools for interpolation, differentiation, and integration used in a Numerical Computing laboratory course. It provides implementations of classical algorithms (forward/central/backward differences, Newton and Lagrange interpolation, Richardson extrapolation, trapezoidal integration, and Simpson's 1/3 and 3/8 rules) plus analyzer code to run experiments and record results.

**Architecture**
- **`src/`**: core implementation.
  - **`src/interpolation/`**: interpolation algorithms and helpers. See [src/interpolation/divided_differences.py](src/interpolation/divided_differences.py#L1-L200) for the Newton divided-differences implementation.
  - **`src/differentiation/`**: finite-difference schemes and Richardson extrapolation.
  - **`src/analyzer/`**: scripts that run experiments, collect outputs, and write the CSV/text files in `output/`.
  - **`src/utils/io_utils.py`**: I/O helpers for reading inputs and writing analysis results.
- **`input/`**: example input parameter files used by the analysis scripts.
- **`output/`**: generated tables, CSVs, and human-readable analysis.

**High-level functionality**
- Interpolate a function given sample points using Newton's divided-differences or Lagrange interpolation.
- Compute numerical derivatives (forward, backward, central) and improve accuracy with Richardson extrapolation.
- Approximate definite integrals with the composite trapezoidal rule and composite Simpson's 1/3 and 3/8 rules.
- Produce difference tables, coefficient lists, and human-readable reports.

**Module responsibilities**
- `src/interpolation/divided_differences.py`: Builds the divided-difference table and evaluates Newton polynomials efficiently.
- `src/interpolation/newton_interpolation.py`: (wraps the divided-differences object to offer a public interpolator API).
- `src/interpolation/lagrange_interpolation.py`: Lagrange form (useful for theory and small-scale interpolation).
- `src/differentiation/numerical_differentiation.py`: Provides unified interfaces for derivative approximations.
- `src/differentiation/richardson_extrapolation.py`: Improves a base finite-difference estimate using extrapolation.
- `src/integration/`: Numerical integration package with reusable data sources, a shared abstract base, the composite trapezoidal rule, and both composite Simpson rules.
- `src/analyzer/*`: Use these to reproduce tables in `output/` for lab reports.

**Numerical integration**

The integration package exposes the composite trapezoidal rule, composite Simpson's 1/3 rule, and composite Simpson's 3/8 rule through `src/integration/__init__.py`. Every method accepts the same `DataSource` abstraction and inherits the common interface from `IntegrationMethod`.

**Object-oriented structure and workflow**

- `DataSource` defines the read-only `x_values` and `y_values` interface. `FunctionDataSource(func, a, b, n)` samples a callable at `n` evenly spaced points; `ArrayDataSource(x_values, y_values)` wraps measured or precomputed values and checks that both arrays have matching lengths and at least two points.
- `IntegrationMethod` validates that its input is a `DataSource`, declares `integrate()` and `method_name`, checks uniform spacing for Simpson's rules, and supplies `result_row(exact=None)` for common result reporting.
- Each concrete rule inherits from `IntegrationMethod` and implements only its numerical formula. A caller can therefore change the rule without changing how the data is supplied or how its result is reported.
- `result_row()` returns an `IntegrationResultRow` with the method name, interval bounds, number of points, and approximation. Passing `exact` also fills in the exact value and absolute error `abs(exact - result)`.

The calculation flow is: create a data source, instantiate a rule with it, call `integrate()` for a float or `result_row(exact=...)` for structured output. Function samples are generated when the `FunctionDataSource` is constructed; array data is not resampled.

**Composite Simpson formulas**

Let `m` be the number of subintervals and `h` the uniform signed spacing between samples. The number of supplied points is `m + 1`. Simpson's rules require equal spacing; they support reversed bounds because the signed spacing is retained.

For Simpson's 1/3 rule, `m` must be even and at least two:

$$\int_a^b f(x)\,dx \approx \frac{h}{3}\left[y_0+y_m+4\sum_{\substack{i=1\\i\text{ odd}}}^{m-1}y_i+2\sum_{\substack{i=2\\i\text{ even}}}^{m-2}y_i\right]$$

`CompositeSimpsonOneThird` applies weights `1, 4, 2, 4, ..., 2, 4, 1` to the sampled values.

For Simpson's 3/8 rule, `m` must be a multiple of three and at least three:

$$\int_a^b f(x)\,dx \approx \frac{3h}{8}\left[y_0+y_m+3\sum_{\substack{1\leq i<m\\i\not\equiv 0\pmod 3}}y_i+2\sum_{\substack{3\leq i<m\\i\equiv 0\pmod 3}}y_i\right]$$

`CompositeSimpsonThreeEighths` applies weight `2` at interior indices divisible by three and weight `3` at the other interior indices. Both rules are exact in exact arithmetic for polynomials up to degree three when their sampling requirements are met.

If a Simpson rule receives too few points, an invalid subinterval count, or nonuniform x-values, `integrate()` raises `ValueError` instead of silently applying an invalid formula. The trapezoidal rule does not require uniform spacing.

**Run integration from the main menu**

From the repository root, start the application and choose `4. Numerical Integration`. The program reads all configured step sizes and automatically compares the trapezoidal, Simpson's 1/3, and Simpson's 3/8 rules on each grid. The provided example uses `h=0.6666666666667, 0.3333333333333, 0.2222222222222` over `[0, 4]`; each step gives a subinterval count divisible by six, so all methods use the same grid at each `h`.

```bash
python main.py
```

For a built-in function, the config uses these keys:

```text
source=function
function=poly    # exp, sin, cos, or poly (x^3 - 2x + 1)
a=0
b=4
h=0.6666666666667, 0.3333333333333, 0.2222222222222
exact=52         # optional; enables absolute-error reporting
```

For measured samples instead, set `source=array` and provide equally spaced values with a number of subintervals divisible by six:

```text
source=array
x=0, 1, 2, 3, 4, 5, 6
y=1, 2, 5, 10, 17, 26, 37
# exact=enter-a-number-here-if-known
```

The `exact` line is optional for either input type. Lines beginning with `#` and blank lines are ignored. Results for all rules and step sizes are compared in `output/integration_comparison_table.txt`, `output/integration_comparison_results.csv`, `output/integration_comparison_analysis.txt`, and `output/integration_comparison_graph.png`. The graph plots each method's approximate integral against `h` and includes the exact integral as a reference when provided.

The Simpson methods validate that x-values are equally spaced and that the subinterval count is valid for the selected rule; the trapezoidal method accepts nonuniform x-values.

**Detailed line-by-line explanation (example):**
The following section explains `src/interpolation/divided_differences.py` line-by-line so you can see why each part exists and how it supports the overall goals.

File: [src/interpolation/divided_differences.py](src/interpolation/divided_differences.py#L1-L200)

- Top docstring and imports:
  - The module docstring documents purpose: build divided-difference table and support Newton interpolation. Clear documentation is critical for maintainability.
  - `from typing import List` brings type annotations; they make interfaces explicit and help static analysis.
  - `from .interpolation_base import InterpolationBase` imports a shared base class that centralizes validation and shared data (`_x_data`, `_y_data`, `_n`). Using a base class avoids duplication and makes behavior consistent across interpolators.

- `class DividedDifferences(InterpolationBase):`
  - Inherits `InterpolationBase` to get input validation and canonical properties (like `_n`). Inheritance clarifies this class is an interpolator type.

- `def __init__(self, x_data: List[float], y_data: List[float]) -> None:`
  - Calls `super().__init__(x_data, y_data)` to reuse validation/assignment performed by the base class. This ensures `x` and `y` arrays are consistent in length and type.
  - Stores the computed divided-difference table in `self._dd_table` by calling `_build_table()` so the heavy computation happens once on construction.

- `_build_table(self) -> List[List[float]]:`
  - `n = self._n` captures the number of points. Local variables reduce attribute lookups and clarify intent.
  - `table: List[List[float]] = [[0.0] * n for _ in range(n)]` allocates a lower-triangular storage as an n-by-n matrix (unused upper entries remain zero). Using a simple list-of-lists is fast, straightforward, and memory is trivial for lab-sized inputs.
  - First loop copies the 0th-order differences (`y` values) into `table[i][0]` — these are the base of the recurrence.
  - Double loop: for each order `j` from 1..n-1 and for rows `i` starting at `j`, compute the divided-difference recurrence
    (table[i][j] - table[i-1][j-1]) / (x[i] - x[i-j]). This implements the formula f[x_i,...,x_{i-j}] = (f[x_i,...,x_{i-j+1}] - f[x_{i-1},...,x_{i-j}])/(x_i - x_{i-j}).
  - Returning `table` keeps the implementation pure: input → computed table.

- `@property def dd_table(self) -> List[List[float]]:`
  - Returns a shallow copy (`[row[:] for row in self._dd_table]`) to prevent callers from mutating internal state. Defensive copies are important for robust APIs.

- `@property def coefficients(self) -> List[float]:`
  - Returns the diagonal entries `self._dd_table[i][i]`, which are the Newton polynomial coefficients `f[x_0,...,x_i]`. Having a `coefficients` accessor makes it easy to evaluate or export the polynomial.

- `def interpolate(self, x: float) -> float:`
  - Evaluates the Newton polynomial using a Horner-like nested multiplication. This is O(n) per evaluation versus a naive O(n^2) expansion.
  - Implementation starts from the highest-order coefficient and iteratively multiplies by `(x - x_k)` and adds the next coefficient.
  - This form is numerically stable (relatively) and efficient; prefer it for repeated evaluations.

- `def format_table(self) -> str:`
  - Creates a human-readable table showing `x_i` and the divided-differences in each column. Useful for lab reports and debugging.

Why this design is chosen (rationale & alternatives):
- The divided-differences table is an industry-standard approach to build Newton polynomials. It separates coefficient computation (one-time O(n^2)) from evaluation (O(n)).
- Storage: a dense n×n list-of-lists is simple; for very large problems you might use sparse or triangular storage to save half the memory, but for educational/numerical-lab sizes this is unnecessary complexity.
- Alternatives:
  - Lagrange form: simpler to write but O(n^2) per evaluation and less convenient for incremental updates.
  - Barycentric Lagrange: superior for numerical stability and repeated evaluation of Lagrange polynomials; consider `barycentric` when interpolating with polynomials based on arbitrary nodes for production code.

Coding best practices applied here
- Input validation centralised in `InterpolationBase` to avoid repeating checks.
- Type annotations for clarity and static checks.
- Small helper methods and properties (`coefficients`, `dd_table`, `format_table`) provide a clean API.
- Defensive copying in accessors to preserve object invariants.

Suggested next steps and learning path
- For class projects: run experiments in `src/analyzer/` with the provided `input/` files and verify `output/` matches expected tables.
- For improved numerical robustness: implement barycentric interpolation and compare against Newton/Lagrange for equispaced and Chebyshev nodes.
- For performance and larger data: consider vectorized NumPy implementations and triangular storage to reduce memory and increase speed.
- Study numerical stability: run experiments plotting error vs. degree for equispaced vs Chebyshev nodes; add notebooks in `Jupyter_notebook-files/` to visualize.

How to run the lab scripts (example)
Run analyzers from the repository root; typical usage (adjust per script):

```bash
python main.py
# or run a specific analyzer, e.g.:
python src/analyzer/interpolation_analyzer.py input/input_interpolation.txt
```

If you want, I can:
- add line-by-line annotated explanations for additional modules (e.g. differentiation files),
- convert key loops to NumPy for speed, or
- add a short `README`-level quickstart script.

---
Generated by assistant to explain architecture, functionality, and coding rationale for this lab repository.

**Quickstart**

Install required packages (use a virtual environment):

```bash
python -m pip install -r requirements.txt
```

Run the non-interactive quickstart to produce output files and plots:

```bash
python quickstart.py
```

Notes on NumPy vectorization
- Where possible the analyzers will attempt to use NumPy to evaluate many
  points or step-sizes at once for speed. This requires that the test
  functions accept NumPy array inputs (i.e., use NumPy ufuncs). If a function
  is a plain Python/math function the code falls back to the scalar
  implementation for full compatibility.

**Differentiation modules — Line-by-line (annotated)**

File: [src/differentiation/numerical_differentiation.py](src/differentiation/numerical_differentiation.py#L1-L200)

- Abstract base class for differentiation methods. Validates inputs (`func` must be callable and `h>0`) ensuring consistent behavior across concrete implementations.
- `derivative(self, x)` is abstract; subclasses provide specific finite-difference formulas.
- `absolute_error(exact, approximate)` is a shared helper for error computation used by analyzers.

File: [src/differentiation/forward_difference.py](src/differentiation/forward_difference.py#L1-L200)

- Implements the forward-difference formula (f(x+h)-f(x))/h as `derivative`.
- This is first-order accurate (error O(h)); simple and inexpensive but less accurate than central differences for the same h.

File: [src/differentiation/backward_difference.py](src/differentiation/backward_difference.py#L1-L200)

- Mirror of the forward difference: (f(x)-f(x-h))/h. Useful when function values ahead of x are unavailable.

File: [src/differentiation/central_difference.py](src/differentiation/central_difference.py#L1-L200)

- Central difference: (f(x+h)-f(x-h))/(2h). Second-order accurate (error O(h^2)) and usually more accurate than forward/backward for smooth functions.

File: [src/differentiation/richardson_extrapolation.py](src/differentiation/richardson_extrapolation.py#L1-L200)

- Uses two central-difference estimates at `h` and `h/2` and combines them as (4*d_half - d_h)/3 to cancel leading-order error, producing a higher-order result. This is a small-step Richardson extrapolation tailored to central differences.

Why these choices matter (rationale):
- Use forward/backward for one-sided derivative estimates or when data is limited to one side.
- Use central for interior points due to better accuracy per h.
- Use Richardson to boost accuracy when function evaluations are cheap and you can compute at `h/2`.

Practical tips
- Always test with decreasing `h` to observe the error plateau: rounding error will dominate for very small `h` values (look for U-shaped error vs h curves).
- When benchmarking, prefer NumPy-aware functions for speed if evaluating over many points or step sizes.


