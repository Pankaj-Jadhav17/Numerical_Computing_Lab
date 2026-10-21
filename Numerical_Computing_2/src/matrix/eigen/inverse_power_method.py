"""Inverse iteration for the smallest-magnitude eigenpair."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .eigen_base import EigenvalueMethod


class InversePowerMethod(EigenvalueMethod):
    """Find the eigenvalue of smallest magnitude using inverse iteration."""

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, object]:
        matrix = self._matrix
        if np.linalg.matrix_rank(matrix) < self.dimension:
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
            eigenvalue = self._rayleigh_quotient(matrix, vector)
            residual = self._residual_norm(matrix, eigenvalue, vector)
            if residual <= self.tolerance:
                return self._single_result(
                    "InversePowerMethod", float(eigenvalue), vector, iteration, True, residual
                )
            last_eigenvalue = eigenvalue

        residual = (
            self._residual_norm(matrix, last_eigenvalue, vector)
            if last_eigenvalue is not None
            else float("inf")
        )
        return self._single_result(
            "InversePowerMethod",
            float(last_eigenvalue) if last_eigenvalue is not None else None,
            vector,
            self.max_iterations,
            False,
            residual,
        )