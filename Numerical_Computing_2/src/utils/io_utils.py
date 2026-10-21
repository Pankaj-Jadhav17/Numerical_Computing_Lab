"""
Input/Output utilities for the Numerical Computing project.

Functions:
  - read_step_size_input(path) -> (x0, [h1, h2, ...])
  - ensure_output_dir(path="output") -> str

The expected input file format is:
  x0
  count_of_h
  h1 h2 h3 ...

Lines starting with '#' or blank lines are ignored.
"""

import os
import math
from typing import List, Tuple

from integration import ArrayDataSource, DataSource, FunctionDataSource


def read_step_size_input(path: str) -> Tuple[float, List[float]]:
    """Read numeric tokens from an input file and return (x0, h_values)."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    tokens: List[float] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            for tok in parts:
                try:
                    tokens.append(float(tok))
                except ValueError as exc:
                    raise ValueError(f"Invalid number in {path}: '{tok}'") from exc

    if len(tokens) < 2:
        raise ValueError(
            f"{path}: expected at least 2 numeric values (x0, count), got {len(tokens)}"
        )

    x0 = tokens[0]
    count = int(tokens[1])
    h_values = tokens[2 : 2 + count]

    if len(h_values) != count:
        raise ValueError(
            f"{path}: header says {count} step sizes but found {len(h_values)}"
        )

    return x0, h_values


def ensure_output_dir(path: str = "output") -> str:
    """Create the output directory if it doesn't exist and return its path."""
    os.makedirs(path, exist_ok=True)
    return path


def read_integration_input(
    path: str,
) -> tuple[DataSource, float | None, list[float]]:
    """Read function integration bounds and steps, or one explicit sample grid."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    config: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                raise ValueError(f"{path}:{line_number}: expected key=value")
            key, value = (part.strip() for part in line.split("=", 1))
            if not key or not value:
                raise ValueError(f"{path}:{line_number}: key and value are required")
            if key in config:
                raise ValueError(f"{path}:{line_number}: duplicate key '{key}'")
            config[key] = value

    source_type = config.get("source", "function").lower()
    allowed = {
        "function": {"source", "function", "a", "b", "h", "exact"},
        "array": {"source", "x", "y", "exact"},
    }
    if source_type not in allowed:
        raise ValueError("source must be 'function' or 'array'")
    unknown = config.keys() - allowed[source_type]
    if unknown:
        raise ValueError(f"Unknown integration input key(s): {', '.join(sorted(unknown))}")

    def numbers(key: str) -> list[float]:
        try:
            values = [float(value) for value in config[key].replace(",", " ").split()]
        except (KeyError, ValueError) as error:
            raise ValueError(f"'{key}' must contain valid numeric value(s)") from error
        if not values or not all(math.isfinite(value) for value in values):
            raise ValueError(f"'{key}' must contain finite numeric value(s)")
        return values

    exact_values = numbers("exact") if "exact" in config else None
    if exact_values is not None and len(exact_values) != 1:
        raise ValueError("'exact' must contain exactly one numeric value")
    exact = exact_values[0] if exact_values else None

    if source_type == "function":
        required = {"function", "a", "b", "h"}
        missing = required - config.keys()
        if missing:
            raise ValueError(f"Missing integration input key(s): {', '.join(sorted(missing))}")
        from math import cos, exp, sin

        functions = {
            "exp": exp,
            "sin": sin,
            "cos": cos,
            "poly": lambda x: x**3 - 2 * x + 1,
        }
        name = config["function"].lower()
        if name not in functions:
            raise ValueError(f"function must be one of: {', '.join(functions)}")
        bounds = numbers("a") + numbers("b")
        h_values = numbers("h")
        if len(bounds) != 2:
            raise ValueError("'a' and 'b' each need one numeric value")
        if bounds[0] == bounds[1]:
            raise ValueError("'a' and 'b' must be different")
        if any(step <= 0 for step in h_values):
            raise ValueError("all h values must be positive")
        if len(set(h_values)) != len(h_values):
            raise ValueError("h values must not contain duplicates")

        interval_length = abs(bounds[1] - bounds[0])
        subinterval_counts = []
        for step in h_values:
            count = interval_length / step
            rounded_count = round(count)
            if not math.isclose(count, rounded_count, rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"interval length must be an integer multiple of h={step}")
            if rounded_count < 6 or rounded_count % 6:
                raise ValueError(
                    f"h={step} must produce a subinterval count divisible by 6 "
                    "for a direct comparison of all three methods"
                )
            subinterval_counts.append(rounded_count)

        data_source = FunctionDataSource(
            functions[name], bounds[0], bounds[1], subinterval_counts[0] + 1
        )
    else:
        missing = {"x", "y"} - config.keys()
        if missing:
            raise ValueError(f"Missing integration input key(s): {', '.join(sorted(missing))}")
        x_values = numbers("x")
        y_values = numbers("y")
        if len(x_values) < 2:
            raise ValueError("array comparison requires at least two data points")
        if len(x_values) != len(y_values):
            raise ValueError("'x' and 'y' must contain the same number of values")
        spacing = x_values[1] - x_values[0]
        if spacing == 0 or any(
            not math.isclose(
                x_value,
                x_values[0] + index * spacing,
                rel_tol=1e-9,
                abs_tol=1e-12,
            )
            for index, x_value in enumerate(x_values[2:], start=2)
        ):
            raise ValueError("array comparison requires equally spaced x values")
        intervals = len(x_values) - 1
        if intervals < 6 or intervals % 6:
            raise ValueError(
                "array comparison requires a subinterval count divisible by 6"
            )
        data_source = ArrayDataSource(x_values, y_values)
        h_values = [abs(spacing)]

    return data_source, exact, h_values


def read_interpolation_input(
    path: str,
) -> Tuple[List[float], List[float], List[float]]:
    """Read an interpolation input file.

    Expected format:
        n                         # number of data points
        x1 x2 ... xn             # x-data values
        [y1 y2 ... yn]           # y-data values (optional blank line)
        m                         # number of query points
        q1 q2 ... qm             # query x-values

    Returns (x_data, y_data, query_points).
    y_data is empty if no y-values were provided.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    data_lines: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(stripped)

    if len(data_lines) < 3:
        raise ValueError(
            f"{path}: expected at least 3 data lines (n, x-values, m/query), "
            f"got {len(data_lines)}"
        )

    n = int(data_lines[0])
    x_data = [float(v) for v in data_lines[1].split()]

    if len(x_data) != n:
        raise ValueError(
            f"{path}: expected {n} x-values, got {len(x_data)}"
        )

    # Determine whether line 3 is y-data or the query-count
    idx = 2
    y_data: List[float] = []

    tokens_line2 = data_lines[idx].split()
    if len(tokens_line2) == n:
        # This line has n values → treat as y-data
        y_data = [float(v) for v in tokens_line2]
        idx += 1

    # Now read query section
    if idx >= len(data_lines):
        raise ValueError(f"{path}: missing query-point count")

    m = int(data_lines[idx])
    idx += 1

    if idx >= len(data_lines):
        raise ValueError(f"{path}: missing query x-values")

    query_points = [float(v) for v in data_lines[idx].split()]
    if len(query_points) != m:
        raise ValueError(
            f"{path}: expected {m} query points, got {len(query_points)}"
        )

    return x_data, y_data, query_points


def read_eigen_input(path: str) -> dict[str, object]:
    """Read an eigenvalue input file through the matrix package handler."""
    from matrix.input_handler import read_eigen_input as read_matrix_input

    return read_matrix_input(path)

