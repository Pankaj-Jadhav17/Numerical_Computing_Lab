"""Validated dense matrix storage and shared operations."""

from __future__ import annotations

from typing import Sequence

import numpy as np

MatrixLike = Sequence[Sequence[float]]


class Matrix:
    """Base class for algorithms that operate on square numeric matrices."""

    def __init__(self, matrix: MatrixLike) -> None:
        try:
            array = np.asarray(matrix, dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError("Matrix must contain finite numeric values") from error

        if array.ndim != 2:
            raise ValueError("Matrix must be a two-dimensional array")
        if array.shape[0] == 0 or array.shape[1] == 0:
            raise ValueError("Matrix must be non-empty")
        if array.shape[0] != array.shape[1]:
            raise ValueError("Matrix must be square")
        if not np.all(np.isfinite(array)):
            raise ValueError("Matrix entries must be finite numbers")
        self._matrix = array

    @property
    def matrix(self) -> np.ndarray:
        return self._matrix.copy()

    @property
    def dimension(self) -> int:
        return self._matrix.shape[0]

    def _matvec(self, vector: np.ndarray) -> np.ndarray:
        return self._matrix @ vector

    @staticmethod
    def _residual_norm(matrix: np.ndarray, eigenvalue: float, eigenvector: np.ndarray) -> float:
        return float(np.linalg.norm(matrix @ eigenvector - eigenvalue * eigenvector))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n={self.dimension})"