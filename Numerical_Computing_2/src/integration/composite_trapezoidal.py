"""Composite trapezoidal rule implementation."""

from .integration_base import IntegrationMethod
from .integration_result import IntegrationResultRow


class CompositeTrapezoidal(IntegrationMethod):
    """Composite trapezoidal rule for integrating sampled function or array data."""

    @property
    def method_name(self) -> str:
        """Human-readable name of the integration method."""
        return "CompositeTrapezoidal"

    def integrate(self) -> float:
        """Approximate the integral using the composite trapezoidal rule."""
        x_values = self.data_source.x_values
        y_values = self.data_source.y_values

        if len(x_values) < 2:
            raise ValueError(
                "At least two points are required for composite trapezoidal integration"
            )

        total = 0.0
        for index in range(len(x_values) - 1):
            x0 = x_values[index]
            x1 = x_values[index + 1]
            y0 = y_values[index]
            y1 = y_values[index + 1]
            total += 0.5 * (y0 + y1) * (x1 - x0)

        return float(total)

    def result_row(self, exact: float | None = None) -> IntegrationResultRow:
        """Create a result row containing the computed integral and optional error."""
        raw_result = self.integrate()
        x_values = self.data_source.x_values

        if exact is None:
            return IntegrationResultRow(
                method=self.method_name,
                interval_start=x_values[0],
                interval_end=x_values[-1],
                n=len(x_values),
                result=raw_result,
            )

        return IntegrationResultRow(
            method=self.method_name,
            interval_start=x_values[0],
            interval_end=x_values[-1],
            n=len(x_values),
            result=raw_result,
            exact=exact,
            absolute_error=self.absolute_error(exact, raw_result),
        )
