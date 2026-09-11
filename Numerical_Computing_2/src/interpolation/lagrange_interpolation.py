"""
Lagrange interpolation implementation.

Provides a numerically stable barycentric form for O(n) evaluation after
an O(n^2) precomputation of barycentric weights.

Class: LagrangeInterpolation(InterpolationBase)
"""
from typing import List

from .interpolation_base import InterpolationBase


class LagrangeInterpolation(InterpolationBase):
    """Lagrange interpolating polynomial using the barycentric formula.

    The constructor precomputes barycentric weights `w_i = 1/\prod_{j!=i}(x_i - x_j)`.
    Evaluations are then performed in O(n) time per point using the formula
    P(x) = (\sum_i w_i y_i /(x - x_i)) / (\sum_i w_i /(x - x_i)).
    """

    def __init__(self, x_data: List[float], y_data: List[float]) -> None:
        super().__init__(x_data, y_data)
        self._weights = self._compute_barycentric_weights()

    def _compute_barycentric_weights(self) -> List[float]:
        n = self._n
        # w_i = 1 / prod_{j != i} (x_i - x_j)
        w: List[float] = [1.0] * n
        for i in range(n):
            prod = 1.0
            xi = self._x_data[i]
            for j in range(n):
                if j == i:
                    continue
                prod *= (xi - self._x_data[j])

            if prod == 0.0:
                # Shouldn't happen because InterpolationBase checks distinct x
                raise ZeroDivisionError("Zero denominator encountered while computing weights")

            w[i] = 1.0 / prod

        return w

    def interpolate(self, x: float) -> float:
        """Evaluate the Lagrange interpolant at x using barycentric formula."""
        # If x matches a node exactly, return the corresponding y to avoid numerical issues
        for xi, yi in zip(self._x_data, self._y_data):
            if x == xi:
                return yi

        numerator = 0.0
        denominator = 0.0
        for wi, xi, yi in zip(self._weights, self._x_data, self._y_data):
            diff = x - xi
            term = wi / diff
            numerator += term * yi
            denominator += term

        return numerator / denominator

    @property
    def method_name(self) -> str:
        return "Lagrange Interpolation"

    @property
    def polynomial_degree(self) -> int:
        return self._n - 1

    def evaluate_range(self, x_start: float, x_end: float, num_points: int = 100) -> List[float]:
        """Evaluate interpolant at `num_points` equally spaced points.

        Returns list of y-values; x-values can be reconstructed by caller.
        """
        if num_points < 2:
            raise ValueError("num_points must be >= 2")

        step = (x_end - x_start) / (num_points - 1)
        x_values = [x_start + i * step for i in range(num_points)]
        return self.interpolate_many(x_values)
