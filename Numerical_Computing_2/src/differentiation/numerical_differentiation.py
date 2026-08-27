from abc import ABC, abstractmethod


class NumericalDifferentiation(ABC):

    def __init__(self, func, h):
        if not callable(func):
            raise TypeError("func must be a function")

        if h <= 0:
            raise ValueError("h must be greater than 0")

        self.func = func
        self.h = h

    @abstractmethod
    def derivative(self, x):
        pass

    @property
    @abstractmethod
    def method_name(self):
        pass

    @staticmethod
    def absolute_error(exact, approximate):
        return abs(exact - approximate)

    def __repr__(self):
        return f"{self.method_name}(h={self.h})"