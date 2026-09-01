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
