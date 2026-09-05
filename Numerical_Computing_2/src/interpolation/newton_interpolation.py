"""
Newton Interpolation — inherits from DividedDifferences.

Adds:
    - add_point(x, y)      → dynamically extend the interpolating polynomial
    - remove_last_point()   → shrink it back
    - evaluate_with_error() → return value + last-term estimate of error
"""

from typing import List, Tuple

from .divided_differences import DividedDifferences


class NewtonInterpolation(DividedDifferences):
    """Full Newton interpolating polynomial with dynamic point management.

    Inherits the divided-difference table from DividedDifferences and adds
    the ability to grow or shrink the dataset on the fly — useful when
    working with large or streaming datasets.
    """

    def __init__(self, x_data: List[float], y_data: List[float]) -> None:
        super().__init__(x_data, y_data)

    # ── dynamic dataset management ─────────────────────────────────────
    def add_point(self, x_new: float, y_new: float) -> None:
        """Append a new data point and update the DD table in O(n) time."""
        if x_new in self._x_data:
            raise ValueError(f"x = {x_new} already exists in the dataset")

        self._x_data.append(x_new)
        self._y_data.append(y_new)
        self._n += 1

        # Extend the DD table with one new row
        new_row = [0.0] * self._n
        new_row[0] = y_new
        idx = self._n - 1

        for j in range(1, self._n):
            new_row[j] = (
                (new_row[j - 1] - self._dd_table[idx - 1][j - 1])
                / (self._x_data[idx] - self._x_data[idx - j])
            )

        self._dd_table.append(new_row)

    def remove_last_point(self) -> Tuple[float, float]:
        """Remove the last data point and shrink the DD table.

        Returns the removed (x, y) pair.
        Raises ValueError if only one point remains.
        """
        if self._n <= 1:
            raise ValueError("Cannot remove the last remaining point")

        x_removed = self._x_data.pop()
        y_removed = self._y_data.pop()
        self._dd_table.pop()
        self._n -= 1

        return x_removed, y_removed

    # ── evaluation with error estimate ─────────────────────────────────
    def evaluate_with_error(self, x: float) -> Tuple[float, float]:
        """Return (interpolated_value, error_estimate).

        The error estimate is |c_n * product(x - x_i)|, where c_n is the
        highest-order divided-difference coefficient.  This is the next
        term in the Newton polynomial and gives a rough indication of the
        interpolation error.
        """
        value = self.interpolate(x)

        # Error estimate: last-term magnitude
        product = 1.0
        for xi in self._x_data:
            product *= (x - xi)

        # The next coefficient would be unknown, so use the last one
        # as an approximation of the magnitude
        last_coeff = self._dd_table[self._n - 1][self._n - 1]
        error_estimate = abs(last_coeff * product)

        return value, error_estimate

    # ── bulk evaluation helpers ────────────────────────────────────────
    def evaluate_range(
        self,
        x_start: float,
        x_end: float,
        num_points: int = 100,
    ) -> Tuple[List[float], List[float]]:
        """Evaluate the polynomial at *num_points* equally-spaced values
        in [x_start, x_end].  Returns (x_values, y_values)."""
        if num_points < 2:
            raise ValueError("num_points must be >= 2")

        step = (x_end - x_start) / (num_points - 1)
        x_values = [x_start + i * step for i in range(num_points)]
        y_values = self.interpolate_many(x_values)

        return x_values, y_values

    @property
    def method_name(self) -> str:
        return "Newton Interpolation"

    @property
    def polynomial_degree(self) -> int:
        """Degree of the interpolating polynomial."""
        return self._n - 1
