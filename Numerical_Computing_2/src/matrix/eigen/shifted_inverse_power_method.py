"""Shifted inverse iteration for an eigenvalue near a target."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .eigen_base import EigenvalueMethod
from ..matrix import MatrixLike


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
        super().__init__(
            matrix,
            tolerance=tolerance,
            max_iterations=max_iterations,
            initial_vector=initial_vector,
        )
        self._mu = float(mu)

    @property
    def mu(self) -> float:
        return self._mu

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, object]:
        matrix = self._matrix
        shifted = matrix - self._mu * np.eye(self.dimension)
        if np.linalg.matrix_rank(shifted) < self.dimension:
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
            eigenvalue = self._rayleigh_quotient(matrix, vector)
            residual = self._residual_norm(matrix, eigenvalue, vector)
            if residual <= self.tolerance:
                return self._single_result(
                    "ShiftedInversePowerMethod",
                    float(eigenvalue),
                    vector,
                    iteration,
                    True,
                    residual,
                    mu=self._mu,
                )
            last_eigenvalue = eigenvalue

        residual = (
            self._residual_norm(matrix, last_eigenvalue, vector)
            if last_eigenvalue is not None
            else float("inf")
        )
        return self._single_result(
            "ShiftedInversePowerMethod",
            float(last_eigenvalue) if last_eigenvalue is not None else None,
            vector,
            self.max_iterations,
            False,
            residual,
            mu=self._mu,
        )