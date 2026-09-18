"""Result data structure for numerical integration methods."""

from dataclasses import dataclass


@dataclass
class IntegrationResultRow:
    """Container for a single numerical integration evaluation."""

    method: str
    interval_start: float
    interval_end: float
    n: int
    result: float
    exact: float | None = None
    absolute_error: float | None = None
