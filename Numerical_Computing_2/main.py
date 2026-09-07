import math
import os
import sys
from typing import Callable

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(PROJECT_DIR, "src")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from analyzer import DifferentiationAnalyzer, InterpolationAnalyzer, RichardsonAnalyzer
from utils import ensure_output_dir, read_interpolation_input, read_step_size_input

TEST_FUNCTIONS = {
    "1": ("exp", math.exp, math.exp),
    "2": ("sin", math.sin, math.cos),
    "3": ("cos", math.cos, lambda x: -math.sin(x)),
    "4": (
        "poly",
        lambda x: x**3 - 2 * x + 1,
        lambda x: 3 * x**2 - 2,
    ),
}

INTERP_FUNCTIONS = {
    "1": ("exp", math.exp),
    "2": ("sin", math.sin),
    "3": ("cos", math.cos),
    "4": ("poly", lambda x: x**3 - 2 * x + 1),
    "5": ("log1p", lambda x: math.log(1 + abs(x))),
}


def print_line(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def choose_functions() -> list[tuple[str, Callable[[float], float], Callable[[float], float]]]:
    print("\nChoose a function:")
    print("1. f(x) = e^x")
    print("2. f(x) = sin(x)")
    print("3. f(x) = cos(x)")
    print("4. f(x) = x^3 - 2x + 1")
    print("5. Use all functions")

    choice = input("Enter your choice: ").strip()

    if choice == "5":
        return list(TEST_FUNCTIONS.values())

    if choice in TEST_FUNCTIONS:
        return [TEST_FUNCTIONS[choice]]

    print("Invalid choice.")
    return []


def choose_interp_functions() -> list[tuple[str, Callable[[float], float]]]:
    """Choose interpolation functions."""
    print("\nChoose a function to interpolate:")
    print("1. f(x) = e^x")
    print("2. f(x) = sin(x)")
    print("3. f(x) = cos(x)")
    print("4. f(x) = x^3 - 2x + 1")
    print("5. f(x) = log(1 + |x|)")
    print("6. Use all functions")

    choice = input("Enter your choice: ").strip()

    if choice == "6":
        return list(INTERP_FUNCTIONS.values())

    if choice in INTERP_FUNCTIONS:
        return [INTERP_FUNCTIONS[choice]]

    print("Invalid choice.")
    return []


def read_input_file(filename: str):
    path = os.path.join(PROJECT_DIR, "input", filename)

    try:
        return read_step_size_input(path)
    except (FileNotFoundError, ValueError) as error:
        print(f"Error: {error}")
        return None


def save_differentiation_results(analyzer: DifferentiationAnalyzer) -> None:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))

    analyzer.write_table(os.path.join(output_dir, "differentiation_table.txt"))
    analyzer.write_csv(os.path.join(output_dir, "differentiation_results.csv"))
    analyzer.write_analysis(os.path.join(output_dir, "differentiation_analysis.txt"))
    analyzer.plot(os.path.join(output_dir, "differentiation_loglog_plot.png"))


def save_richardson_results(analyzer: RichardsonAnalyzer) -> None:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))

    analyzer.write_table(os.path.join(output_dir, "richardson_table.txt"))
    analyzer.write_csv(os.path.join(output_dir, "richardson_results.csv"))
    analyzer.write_analysis(os.path.join(output_dir, "richardson_analysis.txt"))
    analyzer.plot(os.path.join(output_dir, "richardson_loglog_plot.png"))


def save_interpolation_results(analyzer: InterpolationAnalyzer) -> None:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))

    analyzer.write_table(os.path.join(output_dir, "interpolation_table.txt"))
    analyzer.write_csv(os.path.join(output_dir, "interpolation_results.csv"))
    analyzer.write_analysis(os.path.join(output_dir, "interpolation_analysis.txt"))
    analyzer.plot(os.path.join(output_dir, "interpolation_plot.png"))


def run_differentiation() -> None:
    print_line("NUMERICAL DIFFERENTIATION")

    data = read_input_file("input_differentiation.txt")
    if data is None:
        return

    x0, h_values = data
    print(f"Evaluation point: {x0}")
    print("Step sizes:", [f"{h:.0e}" for h in h_values])

    functions = choose_functions()
    if not functions:
        return

    analyzer = DifferentiationAnalyzer(x0, h_values)
    for name, func, exact_derivative in functions:
        analyzer.add_function(name, func, exact_derivative)

    analyzer.run()
    print("\n" + analyzer.format_table())
    print(analyzer.format_analysis())

    save_differentiation_results(analyzer)
    print("\nResults saved in the output folder.")


def run_richardson() -> None:
    print_line("RICHARDSON EXTRAPOLATION")

    data = read_input_file("input_richardson.txt")
    if data is None:
        return

    x0, h_values = data
    print(f"Evaluation point: {x0}")
    print("Step sizes:", [f"{h:.0e}" for h in h_values])

    functions = choose_functions()
    if not functions:
        return

    analyzer = RichardsonAnalyzer(x0, h_values)
    for name, func, exact_derivative in functions:
        analyzer.add_function(name, func, exact_derivative)

    analyzer.run()
    print("\n" + analyzer.format_table())
    print(analyzer.format_analysis())

    save_richardson_results(analyzer)
    print("\nResults saved in the output folder.")


def run_interpolation() -> None:
    print_line("NEWTON INTERPOLATION (DIVIDED DIFFERENCES)")

    print("\nChoose input data source:")
    print("1. Small dataset  (input_interpolation.txt)")
    print("2. Large dataset  (input_interpolation_large.txt)")
    print("3. Enter data manually")

    source_choice = input("Enter your choice: ").strip()

    if source_choice == "1":
        input_file = "input_interpolation.txt"
    elif source_choice == "2":
        input_file = "input_interpolation_large.txt"
    elif source_choice == "3":
        x_data, y_data, query_points = _read_manual_interpolation_input()
        if not x_data:
            return
        _run_interpolation_with_data(x_data, y_data, query_points)
        return
    else:
        print("Invalid choice.")
        return

    path = os.path.join(PROJECT_DIR, "input", input_file)
    try:
        x_data, y_data, query_points = read_interpolation_input(path)
    except (FileNotFoundError, ValueError) as error:
        print(f"Error: {error}")
        return

    _run_interpolation_with_data(x_data, y_data, query_points)


def _read_manual_interpolation_input():
    """Read interpolation data from the user interactively."""
    try:
        n = int(input("Number of data points: ").strip())
        if n < 1:
            print("Need at least 1 data point.")
            return [], [], []

        print(f"Enter {n} x-values (space-separated):")
        x_data = [float(value) for value in input().split()]
        if len(x_data) != n:
            print(f"Expected {n} values, got {len(x_data)}.")
            return [], [], []

        use_function = input("Compute y-values from a function? (y/n): ").strip().lower()
        y_data = []
        if use_function != "y":
            print(f"Enter {n} y-values (space-separated):")
            y_data = [float(value) for value in input().split()]
            if len(y_data) != n:
                print(f"Expected {n} values, got {len(y_data)}.")
                return [], [], []

        m = int(input("Number of query points: ").strip())
        print(f"Enter {m} query x-values (space-separated):")
        query_points = [float(value) for value in input().split()]
        if len(query_points) != m:
            print(f"Expected {m} values, got {len(query_points)}.")
            return [], [], []

        return x_data, y_data, query_points
    except (ValueError, EOFError) as error:
        print(f"Input error: {error}")
        return [], [], []


def _run_interpolation_with_data(x_data: list, y_data: list, query_points: list) -> None:
    """Run interpolation analysis with provided data."""
    print(f"\nData points: {len(x_data)}")
    print(f"x range: [{min(x_data):.4f}, {max(x_data):.4f}]")
    print(f"Query points: {query_points}")

    functions = choose_interp_functions()
    if not functions:
        return

    analyzer = InterpolationAnalyzer(query_points)
    for name, func in functions:
        if y_data:
            analyzer.add_function(name, func, x_data, y_data)
        else:
            analyzer.add_function(name, func, x_data)

    analyzer.run()
    print("\n" + analyzer.format_table())
    print(analyzer.format_dd_tables())
    print(analyzer.format_analysis())

    save_interpolation_results(analyzer)
    print("\nResults saved in the output folder.")


def main() -> None:
    while True:
        print_line("NUMERICAL COMPUTING PROJECT")
        print("1. Numerical Differentiation")
        print("2. Richardson Extrapolation")
        print("3. Newton Interpolation (Divided Differences)")
        print("0. Exit")

        choice = input("Enter your choice: ").strip()

        if choice == "1":
            run_differentiation()
        elif choice == "2":
            run_richardson()
        elif choice == "3":
            run_interpolation()
        elif choice == "0":
            print("Goodbye.")
            break
        else:
            print("Invalid choice.")


if __name__ == "__main__":
    main()
