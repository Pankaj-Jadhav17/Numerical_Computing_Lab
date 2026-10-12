"""Eigenvalue and eigenvector solvers for symmetric and general matrices."""

from __future__ import annotations

import math
import os
from typing import Any, Sequence

import numpy as np

MatrixLike = Sequence[Sequence[float]]


def _as_matrix(matrix: MatrixLike) -> np.ndarray:
    """Convert an input matrix to a validated float array."""
    try:
        array = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Matrix must contain finite numeric values") from exc

    if array.ndim != 2:
        raise ValueError("Matrix must be a two-dimensional array")
    if array.shape[0] == 0 or array.shape[1] == 0:
        raise ValueError("Matrix must be non-empty")
    if array.shape[0] != array.shape[1]:
        raise ValueError("Matrix must be square")
    if not np.all(np.isfinite(array)):
        raise ValueError("Matrix entries must be finite numbers")
    return array


def _normalize(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        raise ValueError("Vector must be non-zero")
    return vector / norm


def _rayleigh_quotient(matrix: np.ndarray, vector: np.ndarray) -> float:
    numerator = float(vector.T @ (matrix @ vector))
    denominator = float(vector.T @ vector)
    if denominator == 0.0:
        raise ValueError("Vector norm is zero")
    return numerator / denominator


def _residual_norm(matrix: np.ndarray, eigenvalue: float, eigenvector: np.ndarray) -> float:
    residual = matrix @ eigenvector - eigenvalue * eigenvector
    return float(np.linalg.norm(residual))


class EigenvalueMethod:
    """Base class for iterative eigenvalue algorithms."""

    def __init__(
        self,
        matrix: MatrixLike,
        *,
        tolerance: float = 1e-10,
        max_iterations: int = 1000,
        initial_vector: Sequence[float] | None = None,
    ) -> None:
        self._matrix = _as_matrix(matrix)
        if tolerance <= 0.0:
            raise ValueError("tolerance must be positive")
        if max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        self._tolerance = float(tolerance)
        self._max_iterations = int(max_iterations)
        self._initial_vector = initial_vector

    @property
    def matrix(self) -> np.ndarray:
        return self._matrix.copy()

    @property
    def tolerance(self) -> float:
        return self._tolerance

    @property
    def max_iterations(self) -> int:
        return self._max_iterations

    def _prepare_vector(self, initial_vector: Sequence[float] | None = None) -> np.ndarray:
        vector_values = (
            np.asarray(initial_vector, dtype=float)
            if initial_vector is not None
            else np.asarray(self._initial_vector, dtype=float)
            if self._initial_vector is not None
            else np.arange(1, self._matrix.shape[0] + 1, dtype=float)
        )
        if vector_values.ndim != 1:
            raise ValueError("Initial vector must be one-dimensional")
        if vector_values.size != self._matrix.shape[0]:
            raise ValueError("Initial vector dimension must match the matrix size")
        if not np.all(np.isfinite(vector_values)):
            raise ValueError("Initial vector entries must be finite numbers")
        vector = np.asarray(vector_values, dtype=float)
        if np.linalg.norm(vector) == 0.0:
            raise ValueError("Initial vector must be non-zero")
        return _normalize(vector)

    @staticmethod
    def _solve_system(matrix: np.ndarray, rhs: np.ndarray, *, label: str) -> np.ndarray:
        try:
            return np.linalg.solve(matrix, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f"{label} requires a nonsingular system") from exc

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n={self._matrix.shape[0]})"


class PowerMethod(EigenvalueMethod):
    """Find the eigenvalue of largest magnitude and its eigenvector."""

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, Any]:
        matrix = self._matrix
        if np.allclose(matrix, 0.0):
            raise ValueError("Power method requires a non-zero matrix")

        vector = self._prepare_vector(initial_vector)
        last_eigenvalue = None

        for iteration in range(1, self.max_iterations + 1):
            next_vector = matrix @ vector
            norm = np.linalg.norm(next_vector)
            if norm == 0.0:
                raise ValueError("Power iteration stalled because A v = 0")
            vector = next_vector / norm
            eigenvalue = _rayleigh_quotient(matrix, vector)
            residual = _residual_norm(matrix, eigenvalue, vector)
            if residual <= self.tolerance:
                return {
                    "method": "PowerMethod",
                    "matrix": matrix,
                    "eigenvalue": float(eigenvalue),
                    "eigenvector": vector,
                    "iterations": iteration,
                    "converged": True,
                    "residual_norm": residual,
                    "tolerance": self.tolerance,
                    "max_iterations": self.max_iterations,
                }
            last_eigenvalue = eigenvalue

        return {
            "method": "PowerMethod",
            "matrix": matrix,
            "eigenvalue": float(last_eigenvalue) if last_eigenvalue is not None else None,
            "eigenvector": vector,
            "iterations": self.max_iterations,
            "converged": False,
            "residual_norm": _residual_norm(matrix, last_eigenvalue, vector)
            if last_eigenvalue is not None
            else float("inf"),
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
        }


class InversePowerMethod(EigenvalueMethod):
    """Find the eigenvalue of smallest magnitude using inverse iteration."""

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, Any]:
        matrix = self._matrix
        if np.linalg.matrix_rank(matrix) < matrix.shape[0]:
            raise ValueError(
                "Inverse power method requires a nonsingular matrix; "
                "zero eigenvalues make the inverse iteration ill-posed"
            )

        vector = self._prepare_vector(initial_vector)
        last_eigenvalue = None

        for iteration in range(1, self.max_iterations + 1):
            next_vector = self._solve_system(matrix, vector, label="Inverse power iteration")
            norm = np.linalg.norm(next_vector)
            if norm == 0.0:
                raise ValueError("Inverse iteration produced a zero vector")
            vector = next_vector / norm
            eigenvalue = _rayleigh_quotient(matrix, vector)
            residual = _residual_norm(matrix, eigenvalue, vector)
            if residual <= self.tolerance:
                return {
                    "method": "InversePowerMethod",
                    "matrix": matrix,
                    "eigenvalue": float(eigenvalue),
                    "eigenvector": vector,
                    "iterations": iteration,
                    "converged": True,
                    "residual_norm": residual,
                    "tolerance": self.tolerance,
                    "max_iterations": self.max_iterations,
                }
            last_eigenvalue = eigenvalue

        return {
            "method": "InversePowerMethod",
            "matrix": matrix,
            "eigenvalue": float(last_eigenvalue) if last_eigenvalue is not None else None,
            "eigenvector": vector,
            "iterations": self.max_iterations,
            "converged": False,
            "residual_norm": _residual_norm(matrix, last_eigenvalue, vector)
            if last_eigenvalue is not None
            else float("inf"),
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
        }


class ShiftedInversePowerMethod(EigenvalueMethod):
    """Find the eigenvalue nearest to a specified target value mu."""

    def __init__(
        self,
        matrix: MatrixLike,
        *,
        mu: float,
        tolerance: float = 1e-10,
        max_iterations: int = 1000,
        initial_vector: Sequence[float] | None = None,
    ) -> None:
        super().__init__(matrix, tolerance=tolerance, max_iterations=max_iterations, initial_vector=initial_vector)
        self._mu = float(mu)

    @property
    def mu(self) -> float:
        return self._mu

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, Any]:
        matrix = self._matrix
        shifted = matrix - self._mu * np.eye(matrix.shape[0])
        if np.linalg.matrix_rank(shifted) < shifted.shape[0]:
            raise ValueError(
                f"Shifted matrix is singular at mu={self._mu}; choose a different target value"
            )

        vector = self._prepare_vector(initial_vector)
        last_eigenvalue = None

        for iteration in range(1, self.max_iterations + 1):
            next_vector = self._solve_system(shifted, vector, label="Shifted inverse iteration")
            norm = np.linalg.norm(next_vector)
            if norm == 0.0:
                raise ValueError("Shifted inverse iteration produced a zero vector")
            vector = next_vector / norm
            eigenvalue = _rayleigh_quotient(matrix, vector)
            residual = _residual_norm(matrix, eigenvalue, vector)
            if residual <= self.tolerance:
                return {
                    "method": "ShiftedInversePowerMethod",
                    "matrix": matrix,
                    "mu": self._mu,
                    "eigenvalue": float(eigenvalue),
                    "eigenvector": vector,
                    "iterations": iteration,
                    "converged": True,
                    "residual_norm": residual,
                    "tolerance": self.tolerance,
                    "max_iterations": self.max_iterations,
                }
            last_eigenvalue = eigenvalue

        return {
            "method": "ShiftedInversePowerMethod",
            "matrix": matrix,
            "mu": self._mu,
            "eigenvalue": float(last_eigenvalue) if last_eigenvalue is not None else None,
            "eigenvector": vector,
            "iterations": self.max_iterations,
            "converged": False,
            "residual_norm": _residual_norm(matrix, last_eigenvalue, vector)
            if last_eigenvalue is not None
            else float("inf"),
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
        }


class JacobiMethod(EigenvalueMethod):
    """Compute all eigenpairs of a real symmetric matrix via cyclic Jacobi rotations."""

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, Any]:
        del initial_vector
        matrix = self._matrix.copy()
        if not np.allclose(matrix, matrix.T, atol=1e-10, rtol=1e-8):
            raise ValueError("Jacobi method requires a real symmetric matrix")

        n = matrix.shape[0]
        if n == 1:
            eigenvalue = float(matrix[0, 0])
            return {
                "method": "JacobiMethod",
                "matrix": matrix,
                "eigenvalues": [eigenvalue],
                "eigenvectors": [np.array([1.0])],
                "iterations": 0,
                "converged": True,
                "residual_norms": [0.0],
                "tolerance": self.tolerance,
                "max_iterations": self.max_iterations,
            }

        a = matrix.copy()
        v = np.eye(n)
        converged = False
        last_max_off = float("inf")

        for iteration in range(1, self.max_iterations + 1):
            p, q = 0, 1
            max_off = 0.0
            for i in range(n):
                for j in range(i + 1, n):
                    value = abs(a[i, j])
                    if value > max_off:
                        max_off = value
                        p, q = i, j
            if max_off <= self.tolerance:
                converged = True
                last_max_off = max_off
                break
            last_max_off = max_off

            if a[p, q] == 0.0:
                break

            tau = (a[q, q] - a[p, p]) / (2.0 * a[p, q])
            t = math.copysign(1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau)
            c = 1.0 / math.sqrt(1.0 + t * t)
            s = t * c

            app = a[p, p]
            aqq = a[q, q]
            apq = a[p, q]

            a[p, p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
            a[q, q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
            a[p, q] = 0.0
            a[q, p] = 0.0

            for r in range(n):
                if r == p or r == q:
                    continue
                arp = a[r, p]
                arq = a[r, q]
                a[r, p] = c * arp - s * arq
                a[r, q] = s * arp + c * arq
                a[p, r] = a[r, p]
                a[q, r] = a[r, q]

            for r in range(n):
                vrp = v[r, p]
                vrq = v[r, q]
                v[r, p] = c * vrp - s * vrq
                v[r, q] = s * vrp + c * vrq

        eigenvalues = np.diag(a)
        order = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[order]
        eigenvectors = v[:, order]

        residual_norms = []
        for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors.T):
            residual_norms.append(_residual_norm(matrix, float(eigenvalue), eigenvector))

        return {
            "method": "JacobiMethod",
            "matrix": matrix,
            "eigenvalues": [float(value) for value in eigenvalues],
            "eigenvectors": [np.asarray(vector, dtype=float) for vector in eigenvectors.T],
            "iterations": iteration if iteration else 0,
            "converged": converged,
            "residual_norms": residual_norms,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
            "max_off_diagonal": last_max_off,
            "approximation_error": float(np.linalg.norm(matrix - eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T)),
        }


def _matrix_from_string(value: str) -> list[list[float]]:
    rows: list[list[float]] = []
    for row_text in value.split(";"):
        row_text = row_text.strip()
        if not row_text:
            continue
        values = [float(token) for token in row_text.replace(",", " ").split()]
        rows.append(values)
    if not rows:
        raise ValueError("Matrix value is required")
    return rows


def _vector_from_string(value: str) -> list[float]:
    text = value.replace(",", " ")
    return [float(token) for token in text.split()]


def read_eigen_input(path: str) -> dict[str, Any]:
    """Read eigenvalue input files.

    Supported syntax includes either key=value lines or a matrix block:

        method=power
        n=2
        matrix
        4 1
        1 3
        initial_vector=1 0
        tolerance=1e-10
        max_iterations=1000

    or

        method=power
        matrix=4 1; 1 3
        initial_vector=1 0
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            processed = line.split("#", 1)[0].strip()
            if not processed:
                continue
            lines.append(processed)

    config: dict[str, Any] = {}
    index = 0
    while index < len(lines):
        line = lines[index]
        if "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            key = key.lower()
            value = value.strip()
            if key in {"matrix", "initial_vector"}:
                if value:
                    config[key] = value
                else:
                    config[key] = None
            else:
                config[key] = value
            index += 1
            continue

        if line.lower() == "matrix":
            matrix_rows: list[list[float]] = []
            n_value = config.get("n")
            expected_rows = int(n_value) if n_value is not None else None
            index += 1
            while index < len(lines):
                candidate = lines[index]
                if "=" in candidate and candidate.split("=", 1)[0].strip().lower() not in {"matrix"}:
                    break
                processed = candidate.split("#", 1)[0].strip()
                if processed:
                    matrix_rows.append([float(number) for number in processed.replace(",", " ").split()])
                    if expected_rows is not None and len(matrix_rows) >= expected_rows:
                        index += 1
                        break
                index += 1
            config["matrix"] = matrix_rows
            continue

        raise ValueError(f"{path}: unsupported input syntax on line '{line}'")

    if "method" not in config:
        raise ValueError("Eigen input must define a method name")

    method_name = str(config["method"]).lower()
    if method_name == "inverse" or method_name == "inversepower":
        method_name = "inverse"
    elif method_name == "shifted" or method_name == "shiftedinverse":
        method_name = "shifted"
    elif method_name == "jacobi":
        method_name = "jacobi"
    elif method_name not in {"power", "inverse", "shifted", "jacobi"}:
        raise ValueError("method must be one of: power, inverse, shifted, jacobi")
    config["method"] = method_name

    if "matrix" not in config:
        raise ValueError("Eigen input must define a matrix")

    matrix_value = config["matrix"]
    if isinstance(matrix_value, str):
        matrix = _matrix_from_string(matrix_value)
    elif isinstance(matrix_value, list):
        matrix = matrix_value
    else:
        raise ValueError("Matrix value is malformed")

    if not matrix or not all(len(row) == len(matrix[0]) for row in matrix):
        raise ValueError("Matrix rows must all have the same length")
    config["matrix"] = matrix

    if "initial_vector" in config and config["initial_vector"] is not None:
        config["initial_vector"] = _vector_from_string(str(config["initial_vector"]))
    if "tolerance" in config:
        config["tolerance"] = float(config["tolerance"])
    else:
        config["tolerance"] = 1e-10
    if "max_iterations" in config:
        config["max_iterations"] = int(config["max_iterations"])
    else:
        config["max_iterations"] = 1000
    if "mu" in config:
        config["mu"] = float(config["mu"])

    return config


def save_eigen_results(result: dict[str, Any], path: str | None = None) -> str:
    """Save a computed eigenvalue result to a text file."""
    if path is None:
        output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "output")
        os.makedirs(output_dir, exist_ok=True)
        method_name = str(result.get("method", "eigen_solver")).lower()
        path = os.path.join(output_dir, f"{method_name}_results.txt")

    matrix = result.get("matrix")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("Eigenvalue / Eigenvector Analysis\n")
        handle.write(f"Method: {result.get('method', 'unknown')}\n")
        handle.write(f"Converged: {result.get('converged', False)}\n")
        handle.write(f"Iterations: {result.get('iterations', 0)}\n")
        if "mu" in result:
            handle.write(f"Target mu: {result['mu']}\n")
        if matrix is not None:
            handle.write("Input matrix:\n")
            for row in np.asarray(matrix):
                handle.write("  [" + " ".join(f"{value: .12g}" for value in row) + "]\n")
        if "eigenvalue" in result:
            handle.write(f"Dominant eigenvalue: {result['eigenvalue']}\n")
            handle.write(f"Residual norm: {result.get('residual_norm', float('nan'))}\n")
            handle.write("Eigenvector:\n")
            for value in np.asarray(result["eigenvector"], dtype=float):
                handle.write(f"  {value: .12g}\n")
        if "eigenvalues" in result:
            handle.write("Eigenvalues:\n")
            for value in result["eigenvalues"]:
                handle.write(f"  {value: .12g}\n")
            handle.write("Eigenvectors:\n")
            for index, vector in enumerate(result["eigenvectors"]):
                handle.write(f"  v{index + 1}: [" + " ".join(f"{value: .12g}" for value in np.asarray(vector, dtype=float)) + "]\n")
            if "residual_norms" in result:
                handle.write("Residual norms:\n")
                for index, value in enumerate(result["residual_norms"]):
                    handle.write(f"  ||A v{index + 1} - lambda{index + 1} v{index + 1}|| = {value: .12g}\n")
            if "approximation_error" in result:
                handle.write(f"V D V^T approximation error: {result['approximation_error']:.12g}\n")

    return path
