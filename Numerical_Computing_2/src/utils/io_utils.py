import os
from typing import List, Tuple


def read_step_size_input(path: str) -> Tuple[float, List[float]]:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    values = []

    try:
        with open(path, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if line and not line.startswith("#"):
                    values.extend(float(value) for value in line.split())
    except ValueError as error:
        raise ValueError(f"Input file contains an invalid number: {error}") from error

    if len(values) < 2:
        raise ValueError("Input file must contain x0 and the number of h values")

    x0 = values[0]
    count = int(values[1])

    if count <= 0:
        raise ValueError("The number of h values must be greater than 0")

    h_values = values[2:]

    if len(h_values) != count:
        raise ValueError(
            f"Expected {count} h values, but found {len(h_values)}"
        )

    if any(h <= 0 for h in h_values):
        raise ValueError("All h values must be greater than 0")

    return x0, h_values


def ensure_output_dir(path: str = "output") -> str:
    os.makedirs(path, exist_ok=True)
    return path
