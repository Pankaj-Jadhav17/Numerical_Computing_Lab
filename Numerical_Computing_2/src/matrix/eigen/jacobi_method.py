"""Cyclic Jacobi rotations for real symmetric matrices."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from .eigen_base import EigenvalueMethod


class JacobiMethod(EigenvalueMethod):
    """Compute all eigenpairs of a real symmetric matrix via cyclic Jacobi rotations."""

    def solve(self, initial_vector: Sequence[float] | None = None) -> dict[str, object]:
        del initial_vector
        matrix = self._matrix.copy()
        if not np.allclose(matrix, matrix.T, atol=1e-10, rtol=1e-8):
            raise ValueError("Jacobi method requires a real symmetric matrix")

        n = self.dimension
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
        vectors = np.eye(n)
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
            app, aqq, apq = a[p, p], a[q, q], a[p, q]
            a[p, p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
            a[q, q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
            a[p, q] = a[q, p] = 0.0

            for r in range(n):
                if r == p or r == q:
                    continue
                arp, arq = a[r, p], a[r, q]
                a[r, p] = a[p, r] = c * arp - s * arq
                a[r, q] = a[q, r] = s * arp + c * arq

            for r in range(n):
                vrp, vrq = vectors[r, p], vectors[r, q]
                vectors[r, p] = c * vrp - s * vrq
                vectors[r, q] = s * vrp + c * vrq

        eigenvalues = np.diag(a)
        order = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[order]
        eigenvectors = vectors[:, order]
        residual_norms = [
            self._residual_norm(matrix, float(value), vector)
            for value, vector in zip(eigenvalues, eigenvectors.T)
        ]
        return {
            "method": "JacobiMethod",
            "matrix": matrix,
            "eigenvalues": [float(value) for value in eigenvalues],
            "eigenvectors": [np.asarray(vector, dtype=float) for vector in eigenvectors.T],
            "iterations": iteration,
            "converged": converged,
            "residual_norms": residual_norms,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
            "max_off_diagonal": last_max_off,
            "approximation_error": float(
                np.linalg.norm(matrix - eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T)
            ),
        }