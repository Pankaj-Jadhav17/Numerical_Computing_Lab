from abc import ABC, abstractmethod
from math import isclose

from .data_source import DataSource
from .integration_result import IntegrationResultRow


class IntegrationMethod(ABC):
    """Base class providing a common interface for numerical integration algorithms."""

    def __init__(self, data_source: DataSource) -> None:
        if not isinstance(data_source, DataSource):
            raise TypeError("data_source must be a DataSource instance")

        self._data_source = data_source

    @property
    def data_source(self) -> DataSource:
        """Return the data source used by the integration method."""
        return self._data_source

    @abstractmethod
    def integrate(self) -> float:
        """Return the approximate integral value."""

    @property
    @abstractmethod
    def method_name(self) -> str:
        """Human-readable name of the integration method."""

    @staticmethod
    def absolute_error(exact: float, approximate: float) -> float:
        """Compute the absolute error between exact and approximate values."""
        return abs(exact - approximate)

    def _uniform_spacing(self) -> float:
        """Return the common x spacing, rejecting data that is not uniform."""
        x_values = self.data_source.x_values
        if len(x_values) < 2:
            raise ValueError("At least two points are required for integration")

        spacing = x_values[1] - x_values[0]
        if spacing == 0 or any(
            not isclose(
                x_value,
                x_values[0] + index * spacing,
                rel_tol=1e-9,
                abs_tol=1e-12,
            )
            for index, x_value in enumerate(x_values[2:], start=2)
        ):
            raise ValueError("Simpson's rules require equally spaced x values")
        return spacing

    def result_row(self, exact: float | None = None) -> IntegrationResultRow:
        """Create a result row containing the computed integral and optional error."""
        x_values = self.data_source.x_values
        result = self.integrate()
        values = {
            "method": self.method_name,
            "interval_start": x_values[0],
            "interval_end": x_values[-1],
            "n": len(x_values),
            "result": result,
        }
        if exact is not None:
            values["exact"] = exact
            values["absolute_error"] = self.absolute_error(exact, result)
        return IntegrationResultRow(**values)

    def __repr__(self) -> str:
        return f"{self.method_name}(n={len(self.data_source)})"
