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
from typing import List, Tuple


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

