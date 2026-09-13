import math
import numpy as np

from differentiation import ForwardDifference, BackwardDifference, CentralDifference
from analyzer.differentiation_analyzer import DifferentiationAnalyzer


def test_basic_differences():
    f = math.sin
    df = math.cos
    x0 = 0.5
    h = 1e-3

    fd = ForwardDifference(f, h).derivative(x0)
    bd = BackwardDifference(f, h).derivative(x0)
    cd = CentralDifference(f, h).derivative(x0)

    assert abs(fd - df(x0)) < 1e-3
    assert abs(bd - df(x0)) < 1e-3
    assert abs(cd - df(x0)) < 1e-6


def test_analyzer_vectorized_path():
    x0 = 0.3
    h_values = [1e-1, 1e-2, 1e-3, 1e-4]
    analyzer = DifferentiationAnalyzer(x0, h_values)

    # Use numpy-aware function to trigger vectorized path
    analyzer.add_function("sin_np", np.sin, np.cos)
    analyzer.run()

    assert len(analyzer.results) == len(h_values)
    for r in analyzer.results:
        assert r.function == "sin_np"
