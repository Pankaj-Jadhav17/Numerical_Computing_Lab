"""Dominant-eigenpair power iteration."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .eigen_base import EigenvalueMethod


class PowerMethod(EigenvalueMethod):
    """Find the eigenvalue of largest magnitude and its eigenvector."""

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, object]:
        matrix = self._matrix
        if np.allclose(matrix, 0.0):
            raise ValueError("Power method requires a non-zero matrix")

        vector = self._prepare_vector(initial_vector)
        last_eigenvalue = None
        for iteration in range(1, self.max_iterations + 1):
            next_vector = self._matvec(vector)
            norm = np.linalg.norm(next_vector)
            if norm == 0.0:
                raise ValueError("Power iteration stalled because A v = 0")
            vector = next_vector / norm
            eigenvalue = self._rayleigh_quotient(matrix, vector)
            residual = self._residual_norm(matrix, eigenvalue, vector)
            if residual <= self.tolerance:
                return self._single_result(
                    "PowerMethod", float(eigenvalue), vector, iteration, True, residual
                )
            last_eigenvalue = eigenvalue

        residual = (
            self._residual_norm(matrix, last_eigenvalue, vector)
            if last_eigenvalue is not None
            else float("inf")
        )
        return self._single_result(
            "PowerMethod",
            float(last_eigenvalue) if last_eigenvalue is not None else None,
            vector,
            self.max_iterations,
            False,
            residual,
        )