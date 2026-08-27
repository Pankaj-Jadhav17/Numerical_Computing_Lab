from .numerical_differentiation import NumericalDifferentiation


class CentralSecondDerivative(NumericalDifferentiation):
    def derivative(self, x: float) -> float:
        return (
            self.func(x + self.h)
            - 2.0 * self.func(x)
            + self.func(x - self.h)
        ) / (self.h ** 2)

    @property
    def method_name(self) -> str:
        return "Central Difference (Second Derivative)"
