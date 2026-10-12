import math

import numpy as np
import pytest

from eigen import (
    InversePowerMethod,
    JacobiMethod,
    PowerMethod,
    ShiftedInversePowerMethod,
    read_eigen_input,
)


def test_power_method_finds_dominant_eigenpair():
    matrix = np.array([[4.0, 1.0], [1.0, 3.0]])
    method = PowerMethod(matrix, tolerance=1e-10, max_iterations=1000)

    result = method.solve()

    assert result["converged"] is True
    assert result["eigenvalue"] == pytest.approx(4.618033988749895, rel=1e-8)
    assert result["residual_norm"] < 1e-8
    assert np.linalg.norm(result["eigenvector"]) > 0.0


def test_inverse_power_method_finds_smallest_magnitude_eigenpair():
    matrix = np.array([[2.0, 0.0], [0.0, 5.0]])
    method = InversePowerMethod(matrix, tolerance=1e-10, max_iterations=1000)

    result = method.solve()

    assert result["converged"] is True
    assert result["eigenvalue"] == pytest.approx(2.0, rel=1e-8)
    assert result["residual_norm"] < 1e-8


def test_shifted_inverse_power_method_finds_target_eigenvalue():
    matrix = np.array([[2.0, 0.0], [0.0, 5.0]])
    method = ShiftedInversePowerMethod(matrix, mu=4.0, tolerance=1e-10, max_iterations=1000)

    result = method.solve()

    assert result["converged"] is True
    assert result["eigenvalue"] == pytest.approx(5.0, rel=1e-8)
    assert result["residual_norm"] < 1e-8


def test_jacobi_method_recovers_symmetric_eigenpairs():
    matrix = np.array([
        [4.0, 1.0, 1.0],
        [1.0, 3.0, 0.0],
        [1.0, 0.0, 2.0],
    ])
    method = JacobiMethod(matrix, tolerance=1e-12, max_iterations=1000)

    result = method.solve()
    expected = np.linalg.eigvalsh(matrix)

    assert result["converged"] is True
    assert result["eigenvalues"] == pytest.approx(expected, rel=1e-8)
    assert len(result["eigenvectors"]) == 3
    for eigenvalue, eigenvector in zip(result["eigenvalues"], result["eigenvectors"]):
        residual = np.linalg.norm(matrix @ eigenvector - eigenvalue * eigenvector)
        assert residual < 1e-8


def test_invalid_matrix_shapes_and_singular_cases():
    with pytest.raises(ValueError, match="Matrix must|square"):
        PowerMethod([[1.0, 2.0], [3.0]])

    with pytest.raises(ValueError, match="non-zero"):
        PowerMethod(np.zeros((2, 2))).solve()

    with pytest.raises(ValueError, match="nonsingular"):
        InversePowerMethod(np.zeros((2, 2))).solve()

    with pytest.raises(ValueError, match="symmetric"):
        JacobiMethod([[1.0, 2.0], [3.0, 4.0]]).solve()


def test_read_eigen_input_parses_matrix_and_parameters(tmp_path):
    input_file = tmp_path / "eigen.txt"
    input_file.write_text(
        "method=power\n"
        "n=2\n"
        "matrix\n"
        "4 1\n"
        "1 3\n"
        "initial_vector=1 0\n"
        "tolerance=1e-10\n"
        "max_iterations=100\n",
        encoding="utf-8",
    )

    config = read_eigen_input(str(input_file))

    assert config["method"] == "power"
    assert np.asarray(config["matrix"]) == pytest.approx(np.array([[4.0, 1.0], [1.0, 3.0]]))
    assert np.asarray(config["initial_vector"]) == pytest.approx(np.array([1.0, 0.0]))
    assert config["tolerance"] == pytest.approx(1e-10)
    assert config["max_iterations"] == 100
