from .numerical_differentiation import NumericalDifferentiation


class ForwardDifference(NumericalDifferentiation):
    def derivative(self, x: float) -> float:
        return (self.func(x + self.h) - self.func(x)) / self.h

    @property
    def method_name(self) -> str:
        return "Forward Difference"
