from .numerical_differentiation import NumericalDifferentiation


class BackwardDifference(NumericalDifferentiation):
    def derivative(self, x: float) -> float:
        return (self.func(x) - self.func(x - self.h)) / self.h

    @property
    def method_name(self) -> str:
        return "Backward Difference"
