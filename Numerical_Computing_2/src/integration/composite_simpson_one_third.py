"""Composite Simpson's 1/3 rule implementation."""

from .integration_base import IntegrationMethod


class CompositeSimpsonOneThird(IntegrationMethod):
    """Composite Simpson's 1/3 rule for equally spaced sampled data."""

    @property
    def method_name(self) -> str:
        """Human-readable name of the integration method."""
        return "CompositeSimpsonOneThird"

    def integrate(self) -> float:
        """Approximate the integral using an even number of subintervals."""
        x_values = self.data_source.x_values
        y_values = self.data_source.y_values
        intervals = len(x_values) - 1
        if intervals < 2 or intervals % 2:
            raise ValueError("Simpson's 1/3 rule requires an even number of subintervals")

        spacing = self._uniform_spacing()
        weighted_sum = y_values[0] + y_values[-1]
        weighted_sum += 4 * sum(y_values[index] for index in range(1, intervals, 2))
        weighted_sum += 2 * sum(y_values[index] for index in range(2, intervals, 2))
        return float(spacing * weighted_sum / 3)