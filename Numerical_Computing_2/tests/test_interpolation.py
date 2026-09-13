import math
import random

from interpolation import NewtonInterpolation


def test_newton_interpolation_exact_points():
    # Simple polynomial data: f(x)=x^3-2x+1
    x = [i for i in range(-5, 6)]
    y = [xi ** 3 - 2 * xi + 1 for xi in x]
    interp = NewtonInterpolation(x, y)

    # Interpolant should reproduce the original y-values exactly
    for xi, yi in zip(x, y):
        assert abs(interp.interpolate(xi) - yi) < 1e-12


def test_newton_interpolation_large_dataset():
    # Large dataset but smooth function (sin); check reasonable runtime
    n = 1000
    xs = [i * 0.01 for i in range(n)]
    ys = [math.sin(xi) for xi in xs]
    interp = NewtonInterpolation(xs, ys)

    # Random sample checks
    for _ in range(10):
        q = random.uniform(xs[0], xs[-1])
        approx = interp.interpolate(q)
        exact = math.sin(q)
        assert abs(approx - exact) < 1e-2  # loose tolerance for large-degree polynomial
