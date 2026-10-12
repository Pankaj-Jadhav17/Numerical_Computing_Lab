"""Composite Simpson's 3/8 rule implementation."""

from .integration_base import IntegrationMethod


class CompositeSimpsonThreeEighths(IntegrationMethod):
    """Composite Simpson's 3/8 rule for equally spaced sampled data."""

    @property
    def method_name(self) -> str:
        """Human-readable name of the integration method."""
        return "CompositeSimpsonThreeEighths"

    def integrate(self) -> float:
        """Approximate the integral using a multiple of three subintervals."""
        x_values = self.data_source.x_values
        y_values = self.data_source.y_values
        intervals = len(x_values) - 1
        if intervals < 3 or intervals % 3:
            raise ValueError("Simpson's 3/8 rule requires a multiple of three subintervals")

        spacing = self._uniform_spacing()
        weighted_sum = y_values[0] + y_values[-1]
        weighted_sum += 3 * sum(
            y_values[index]
            for index in range(1, intervals)
            if index % 3
        )
        weighted_sum += 2 * sum(
            y_values[index]
            for index in range(3, intervals, 3)
        )
        return float(3 * spacing * weighted_sum / 8)