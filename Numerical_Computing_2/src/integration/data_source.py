"""Data source abstractions for numerical integration."""

from abc import ABC, abstractmethod


class DataSource(ABC):
    """Abstract interface for function-based or array-based integration data."""

    @property
    @abstractmethod
    def x_values(self) -> list[float]:
        """Return the x coordinates for the data source."""

    @property
    @abstractmethod
    def y_values(self) -> list[float]:
        """Return the y values for the data source."""

    def __len__(self) -> int:
        return len(self.x_values)

    def __iter__(self):
        return iter(zip(self.x_values, self.y_values))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n={len(self)})"


class FunctionDataSource(DataSource):
    """Generate integration data from a callable function on a finite interval."""

    def __init__(self, func, a: float, b: float, n: int) -> None:
        if not callable(func):
            raise TypeError("func must be a function")

        if n < 2:
            raise ValueError("n must be at least 2")

        if a == b:
            raise ValueError("a and b must be different")

        self._func = func
        self._a = float(a)
        self._b = float(b)
        self._n = int(n)
        self._x_values = [
            self._a + (self._b - self._a) * index / (self._n - 1)
            for index in range(self._n)
        ]
        self._y_values = [float(self._func(x)) for x in self._x_values]

    @property
    def func(self):
        """Return the underlying callable function."""
        return self._func

    @property
    def a(self) -> float:
        """Return the left endpoint of the interval."""
        return self._a

    @property
    def b(self) -> float:
        """Return the right endpoint of the interval."""
        return self._b

    @property
    def n(self) -> int:
        """Return the number of sampling points."""
        return self._n

    @property
    def x_values(self) -> list[float]:
        """Return a copy of the x values."""
        return list(self._x_values)

    @property
    def y_values(self) -> list[float]:
        """Return a copy of the y values."""
        return list(self._y_values)

    def evaluate(self, x: float) -> float:
        """Evaluate the underlying function at x."""
        return float(self._func(x))


class ArrayDataSource(DataSource):
    """Represent integration data supplied as explicit arrays of x and y values."""

    def __init__(self, x_values: list[float], y_values: list[float]) -> None:
        if len(x_values) != len(y_values):
            raise ValueError(
                f"x_values and y_values must have the same length, "
                f"got {len(x_values)} and {len(y_values)}"
            )

        if len(x_values) < 2:
            raise ValueError("At least two data points are required")

        self._x_values = [float(value) for value in x_values]
        self._y_values = [float(value) for value in y_values]

    @property
    def x_values(self) -> list[float]:
        """Return a copy of the x values."""
        return list(self._x_values)

    @property
    def y_values(self) -> list[float]:
        """Return a copy of the y values."""
        return list(self._y_values)
