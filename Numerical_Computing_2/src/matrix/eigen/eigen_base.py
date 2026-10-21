"""Shared settings and numerical helpers for eigenvalue methods."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ..matrix import Matrix, MatrixLike


class EigenvalueMethod(Matrix):
    """Base class for iterative eigenvalue algorithms."""

    def __init__(
        self,
        matrix: MatrixLike,
        *,
        tolerance: float = 1e-10,
        max_iterations: int = 1000,
        initial_vector: Sequence[float] | None = None,
    ) -> None:
        super().__init__(matrix)
        if tolerance <= 0.0:
            raise ValueError("tolerance must be positive")
        if max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        self._tolerance = float(tolerance)
        self._max_iterations = int(max_iterations)
        self._initial_vector = initial_vector

    @property
    def tolerance(self) -> float:
        return self._tolerance

    @property
    def max_iterations(self) -> int:
        return self._max_iterations

    @staticmethod
    def _normalize(vector: np.ndarray) -> np.ndarray:
        norm = float(np.linalg.norm(vector))
        if norm == 0.0:
            raise ValueError("Vector must be non-zero")
        return vector / norm

    def _prepare_vector(self, initial_vector: Sequence[float] | None = None) -> np.ndarray:
        vector_values = (
            np.asarray(initial_vector, dtype=float)
            if initial_vector is not None
            else np.asarray(self._initial_vector, dtype=float)
            if self._initial_vector is not None
            else np.arange(1, self.dimension + 1, dtype=float)
        )
        if vector_values.ndim != 1:
            raise ValueError("Initial vector must be one-dimensional")
        if vector_values.size != self.dimension:
            raise ValueError("Initial vector dimension must match the matrix size")
        if not np.all(np.isfinite(vector_values)):
            raise ValueError("Initial vector entries must be finite numbers")
        vector = np.asarray(vector_values, dtype=float)
        if np.linalg.norm(vector) == 0.0:
            raise ValueError("Initial vector must be non-zero")
        return self._normalize(vector)

    @staticmethod
    def _rayleigh_quotient(matrix: np.ndarray, vector: np.ndarray) -> float:
        numerator = float(vector.T @ (matrix @ vector))
        denominator = float(vector.T @ vector)
        if denominator == 0.0:
            raise ValueError("Vector norm is zero")
        return numerator / denominator

    @staticmethod
    def _solve_system(matrix: np.ndarray, rhs: np.ndarray, *, label: str) -> np.ndarray:
        try:
            return np.linalg.solve(matrix, rhs)
        except np.linalg.LinAlgError as error:
            raise ValueError(f"{label} requires a nonsingular system") from error

    def _single_result(
        self,
        method: str,
        eigenvalue: float | None,
        eigenvector: np.ndarray,
        iterations: int,
        converged: bool,
        residual_norm: float,
        *,
        mu: float | None = None,
    ) -> dict[str, object]:
        result: dict[str, object] = {
            "method": method,
            "matrix": self._matrix,
            "eigenvalue": eigenvalue,
            "eigenvector": eigenvector,
            "iterations": iterations,
            "converged": converged,
            "residual_norm": residual_norm,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
        }
        if mu is not None:
            result["mu"] = mu
        return result