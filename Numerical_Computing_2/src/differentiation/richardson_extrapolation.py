from .central_difference import CentralDifference
from .numerical_differentiation import NumericalDifferentiation


class RichardsonExtrapolation(NumericalDifferentiation):
    def central_difference(self, x: float, h: float) -> float:
        return CentralDifference(self.func, h).derivative(x)

    def derivative(self, x: float) -> float:
        d_h = self.central_difference(x, self.h)
        d_half = self.central_difference(x, self.h / 2.0)
        return (4.0 * d_half - d_h) / 3.0

    @property
    def method_name(self) -> str:
        return "Richardson Extrapolation"
