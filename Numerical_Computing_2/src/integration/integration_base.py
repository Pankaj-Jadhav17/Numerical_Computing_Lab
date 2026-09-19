from abc import ABC, abstractmethod

from .data_source import DataSource


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

    def __repr__(self) -> str:
        return f"{self.method_name}(n={len(self.data_source)})"
