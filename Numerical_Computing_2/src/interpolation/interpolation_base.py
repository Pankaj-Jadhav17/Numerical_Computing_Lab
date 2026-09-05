"""
Abstract base class for all interpolation methods.

Subclasses must implement:
    - interpolate(x)  → float
    - method_name     → str (property)
"""

from abc import ABC, abstractmethod
from typing import List


class InterpolationBase(ABC):
    """Base class providing a common interface for interpolation algorithms."""

    def __init__(self, x_data: List[float], y_data: List[float]) -> None:
        if len(x_data) != len(y_data):
            raise ValueError(
                f"x_data and y_data must have the same length, "
                f"got {len(x_data)} and {len(y_data)}"
            )

        if len(x_data) < 1:
            raise ValueError("At least one data point is required")

        # Check for duplicate x values
        if len(set(x_data)) != len(x_data):
            raise ValueError("x_data must contain distinct values")

        self._x_data = list(x_data)
        self._y_data = list(y_data)
        self._n = len(x_data)

    # ── public properties ──────────────────────────────────────────────
    @property
    def x_data(self) -> List[float]:
        """Return a copy of the x data points."""
        return list(self._x_data)

    @property
    def y_data(self) -> List[float]:
        """Return a copy of the y data points."""
        return list(self._y_data)

    @property
    def n(self) -> int:
        """Number of data points."""
        return self._n

    # ── abstract interface ─────────────────────────────────────────────
    @abstractmethod
    def interpolate(self, x: float) -> float:
        """Return the interpolated value at *x*."""

    @property
    @abstractmethod
    def method_name(self) -> str:
        """Human-readable name of the interpolation method."""

    # ── convenience ────────────────────────────────────────────────────
    def interpolate_many(self, x_values: List[float]) -> List[float]:
        """Interpolate at several x-values."""
        return [self.interpolate(x) for x in x_values]

    @staticmethod
    def absolute_error(exact: float, approximate: float) -> float:
        """Compute the absolute error between *exact* and *approximate*."""
        return abs(exact - approximate)

    def __repr__(self) -> str:
        return f"{self.method_name}(n={self._n})"
