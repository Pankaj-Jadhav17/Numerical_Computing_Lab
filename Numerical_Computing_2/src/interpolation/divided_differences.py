"""
Divided Differences table computation.

Inherits from InterpolationBase and builds the Newton divided-difference
table.  Can be used standalone to inspect the DD table or as a parent
class for NewtonInterpolation.
"""

from typing import List

from .interpolation_base import InterpolationBase


class DividedDifferences(InterpolationBase):
    """Compute and store the divided-difference table for a set of data points.

    The table is a lower-triangular matrix stored as a list of lists:
        table[i][j] = f[x_i, x_{i-1}, ..., x_{i-j}]
    where j = 0 is simply y_i.
    """

    def __init__(self, x_data: List[float], y_data: List[float]) -> None:
        super().__init__(x_data, y_data)
        self._dd_table: List[List[float]] = self._build_table()

    # ── core computation ───────────────────────────────────────────────
    def _build_table(self) -> List[List[float]]:
        """Build the divided-difference table using the standard algorithm.

        Works efficiently for both small and very large datasets (O(n^2)
        time and space).
        """
        n = self._n
        # table[i][j]: j-th order divided difference starting at index i
        table: List[List[float]] = [[0.0] * n for _ in range(n)]

        # 0th-order divided differences are just the y values
        for i in range(n):
            table[i][0] = self._y_data[i]

        # Higher-order divided differences
        for j in range(1, n):
            for i in range(j, n):
                table[i][j] = (
                    (table[i][j - 1] - table[i - 1][j - 1])
                    / (self._x_data[i] - self._x_data[i - j])
                )

        return table

    # ── public accessors ───────────────────────────────────────────────
    @property
    def dd_table(self) -> List[List[float]]:
        """Return a copy of the divided-difference table."""
        return [row[:] for row in self._dd_table]

    @property
    def coefficients(self) -> List[float]:
        """Return the diagonal of the DD table (Newton polynomial coefficients).

        coefficients[k] = f[x_0, x_1, ..., x_k]
        """
        return [self._dd_table[i][i] for i in range(self._n)]

    # ── InterpolationBase interface ────────────────────────────────────
    def interpolate(self, x: float) -> float:
        """Evaluate the Newton interpolating polynomial at *x* using
        Horner-like nested multiplication (O(n) per evaluation)."""
        coeffs = self.coefficients
        n = self._n

        # Horner's form: start from the highest-order term
        result = coeffs[n - 1]
        for k in range(n - 2, -1, -1):
            result = result * (x - self._x_data[k]) + coeffs[k]

        return result

    @property
    def method_name(self) -> str:
        return "Divided Differences"

    # ── pretty-printing ────────────────────────────────────────────────
    def format_table(self) -> str:
        """Return a human-readable string of the divided-difference table."""
        lines: List[str] = []
        header_parts = ["  i", f"{'x_i':>10}", f"{'f[.]':>14}"]
        for j in range(1, self._n):
            header_parts.append(f"{'f[' + '.'*(j+1) + ']':>14}")
        lines.append("".join(header_parts))
        lines.append("-" * len(lines[0]))

        for i in range(self._n):
            parts = [f"{i:>3}", f"{self._x_data[i]:>10.4f}"]
            for j in range(i + 1):
                parts.append(f"{self._dd_table[i][j]:>14.6e}")
            lines.append("".join(parts))

        return "\n".join(lines)

    def __str__(self) -> str:
        return self.format_table()
