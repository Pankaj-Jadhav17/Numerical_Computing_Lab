from .numerical_differentiation import NumericalDifferentiation


class CentralDifference(NumericalDifferentiation):
    def derivative(self, x: float) -> float:
        return (self.func(x + self.h) - self.func(x - self.h)) / (2.0 * self.h)

    @property
    def method_name(self) -> str:
        return "Central Difference"
