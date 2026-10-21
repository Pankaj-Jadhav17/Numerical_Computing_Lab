"""Input parsing and result serialization for matrix algorithms."""

from __future__ import annotations

import os
from typing import Any

import numpy as np


def _parse_numeric_list(raw_value: str) -> list[float]:
    values = [
        float(token)
        for token in raw_value.replace(";", " ").replace(",", " ").split()
    ]
    if not values:
        raise ValueError("Expected at least one numeric value")
    return values


def read_eigen_input(path: str) -> dict[str, Any]:
    """Read an eigenvalue input file containing a matrix and method settings."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    lines = []
    with open(path, "r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.split("#", 1)[0].strip()
            if line:
                lines.append(line)

    config: dict[str, Any] = {}
    index = 0
    while index < len(lines):
        line = lines[index]
        if "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            key = key.lower()
            config[key] = value if value else None
            index += 1
            continue
        if line.lower() != "matrix":
            raise ValueError(f"{path}: unexpected input syntax: {line}")

        index += 1
        rows = []
        expected_rows = int(config["n"]) if config.get("n") else None
        while index < len(lines) and "=" not in lines[index]:
            rows.append(_parse_numeric_list(lines[index]))
            index += 1
            if expected_rows is not None and len(rows) >= expected_rows:
                break
        if not rows:
            raise ValueError(f"{path}: matrix rows are missing")
        if expected_rows is not None and len(rows) != expected_rows:
            raise ValueError(
                f"{path}: expected {expected_rows} matrix rows but found {len(rows)}"
            )
        config["matrix"] = rows

    if "method" not in config:
        raise ValueError("Eigen input must define a method")
    method = str(config["method"]).lower()
    aliases = {
        "inversepower": "inverse",
        "inverse_power": "inverse",
        "shiftedinverse": "shifted",
        "shifted_inverse": "shifted",
    }
    method = aliases.get(method, method)
    if method not in {"power", "inverse", "shifted", "jacobi"}:
        raise ValueError("method must be one of: power, inverse, shifted, jacobi")
    config["method"] = method

    if "matrix" not in config:
        raise ValueError("Eigen input must define a matrix")
    matrix_value = config["matrix"]
    if isinstance(matrix_value, str):
        rows = [
            _parse_numeric_list(row)
            for row in matrix_value.split(";")
            if row.strip()
        ]
        config["matrix"] = rows
    elif not isinstance(matrix_value, list):
        raise ValueError("Matrix entries are malformed")

    matrix = config["matrix"]
    if not matrix:
        raise ValueError("Matrix is empty")
    row_lengths = {len(row) for row in matrix}
    if len(row_lengths) != 1:
        raise ValueError("Matrix rows must all have the same length")
    if len(matrix) != len(matrix[0]):
        raise ValueError("Matrix must be square")

    if config.get("initial_vector") is not None:
        config["initial_vector"] = _parse_numeric_list(str(config["initial_vector"]))
    config["tolerance"] = float(config.get("tolerance") or 1e-10)
    config["max_iterations"] = int(config.get("max_iterations") or 1000)
    if config.get("mu") is not None:
        config["mu"] = float(config["mu"])
    return config


def save_eigen_results(result: dict[str, Any], path: str | None = None) -> str:
    """Save a computed eigenvalue result as a human-readable text report."""
    if path is None:
        project_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        output_dir = os.path.join(project_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        method_name = str(result.get("method", "eigen_solver")).lower()
        path = os.path.join(output_dir, f"{method_name}_results.txt")

    matrix = result.get("matrix")
    with open(path, "w", encoding="utf-8") as file:
        file.write("Eigenvalue / Eigenvector Analysis\n")
        file.write(f"Method: {result.get('method', 'unknown')}\n")
        file.write(f"Converged: {result.get('converged', False)}\n")
        file.write(f"Iterations: {result.get('iterations', 0)}\n")
        if "mu" in result:
            file.write(f"Target mu: {result['mu']}\n")
        if matrix is not None:
            file.write("Input matrix:\n")
            for row in np.asarray(matrix):
                file.write("  [" + " ".join(f"{value: .12g}" for value in row) + "]\n")
        if "eigenvalue" in result:
            file.write(f"Dominant eigenvalue: {result['eigenvalue']}\n")
            file.write(f"Residual norm: {result.get('residual_norm', float('nan'))}\n")
            file.write("Eigenvector:\n")
            for value in np.asarray(result["eigenvector"], dtype=float):
                file.write(f"  {value: .12g}\n")
        if "eigenvalues" in result:
            file.write("Eigenvalues:\n")
            for value in result["eigenvalues"]:
                file.write(f"  {value: .12g}\n")
            file.write("Eigenvectors:\n")
            for index, vector in enumerate(result["eigenvectors"]):
                values = " ".join(f"{value: .12g}" for value in np.asarray(vector, dtype=float))
                file.write(f"  v{index + 1}: [{values}]\n")
            if "residual_norms" in result:
                file.write("Residual norms:\n")
                for index, value in enumerate(result["residual_norms"]):
                    file.write(
                        f"  ||A v{index + 1} - lambda{index + 1} v{index + 1}|| = "
                        f"{value: .12g}\n"
                    )
            if "approximation_error" in result:
                file.write(
                    f"V D V^T approximation error: {result['approximation_error']:.12g}\n"
                )
    return path