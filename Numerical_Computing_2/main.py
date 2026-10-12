import csv
import math
import os
import sys
from dataclasses import asdict
from typing import Callable

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(PROJECT_DIR, "src")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from analyzer import DifferentiationAnalyzer, InterpolationAnalyzer, RichardsonAnalyzer
from eigen import (
    InversePowerMethod,
    JacobiMethod,
    PowerMethod,
    ShiftedInversePowerMethod,
    read_eigen_input,
    save_eigen_results,
)
from integration import (
    CompositeSimpsonOneThird,
    CompositeSimpsonThreeEighths,
    CompositeTrapezoidal,
    FunctionDataSource,
)
from utils import (
    ensure_output_dir,
    read_eigen_input as read_eigen_file,
    read_integration_input,
    read_interpolation_input,
    read_step_size_input,
)

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


def save_integration_result(result, data_source) -> tuple[str, str, str, str]:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))
    filename_prefix = {
        "CompositeTrapezoidal": "trapezoidal",
        "CompositeSimpsonOneThird": "simpson_1_3",
        "CompositeSimpsonThreeEighths": "simpson_3_8",
    }.get(result.method)
    if filename_prefix is None:
        raise ValueError(f"Unsupported integration method: {result.method}")

    table_path = os.path.join(output_dir, f"{filename_prefix}_table.txt")
    csv_path = os.path.join(output_dir, f"{filename_prefix}_results.csv")
    analysis_path = os.path.join(output_dir, f"{filename_prefix}_analysis.txt")
    graph_path = os.path.join(output_dir, f"{filename_prefix}_graph.png")

    result_data = asdict(result)
    with open(table_path, "w", encoding="utf-8") as file:
        file.write("Numerical Integration Results\n")
        file.write(f"Method: {result.method}\n")
        file.write(f"Interval: [{result.interval_start}, {result.interval_end}]\n")
        file.write(f"Data points: {result.n}\n")
        file.write("\n")
        file.write(
            f"{'Approximation':<20}{'Exact':<20}{'Absolute Error':<20}\n"
        )
        if result.exact is not None:
            file.write(
                f"{result.result:<20.12g}{result.exact:<20.12g}"
                f"{result.absolute_error:<20.12g}\n"
            )
        else:
            file.write(f"{result.result:<20.12g}{'N/A':<20}{'N/A':<20}\n")

    with open(csv_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(result_data))
        writer.writeheader()
        writer.writerow(result_data)

    with open(analysis_path, "w", encoding="utf-8") as file:
        file.write("Numerical Integration Analysis\n")
        file.write(f"{result.method} used {result.n} data points.\n")
        file.write(
            f"Approximate integral over [{result.interval_start}, "
            f"{result.interval_end}] = {result.result:.12g}.\n"
        )
        if result.absolute_error is not None:
            file.write(f"Absolute error = {result.absolute_error:.12g}.\n")

    save_integration_graph(data_source, result, graph_path)
    return table_path, csv_path, analysis_path, graph_path


def save_integration_graph(data_source, result, path: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    x_values = data_source.x_values
    y_values = data_source.y_values
    if result.method == "CompositeTrapezoidal":
        points_per_panel = 2
        title = "Composite Trapezoidal Rule"
    elif result.method == "CompositeSimpsonOneThird":
        points_per_panel = 3
        title = "Composite Simpson's 1/3 Rule"
    else:
        points_per_panel = 4
        title = "Composite Simpson's 3/8 Rule"

    curve_x = []
    curve_y = []
    for start in range(0, len(x_values) - 1, points_per_panel - 1):
        panel_x = np.asarray(x_values[start : start + points_per_panel])
        panel_y = np.asarray(y_values[start : start + points_per_panel])
        degree = len(panel_x) - 1
        coefficients = np.polyfit(panel_x, panel_y, degree)
        sample_x = np.linspace(float(np.min(panel_x)), float(np.max(panel_x)), 100)
        sample_y = np.polyval(coefficients, sample_x)
        curve_x.extend(sample_x.tolist())
        curve_y.extend(sample_y.tolist())

    curve_points = sorted(zip(curve_x, curve_y))
    curve_x, curve_y = zip(*curve_points)

    figure, axis = plt.subplots(figsize=(8, 5))
    axis.plot(curve_x, curve_y, color="#176B87", label="Rule approximation")
    axis.fill_between(curve_x, curve_y, 0, color="#64CCC5", alpha=0.28)
    axis.scatter(x_values, y_values, color="#D95D39", zorder=3, label="Input samples")
    axis.axhline(0, color="#333333", linewidth=0.8)
    axis.set_title(f"{title}\nIntegral ≈ {result.result:.8g}")
    axis.set_xlabel("x")
    axis.set_ylabel("f(x)")
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def save_integration_comparison(
    rows: list[dict[str, object]], data_source
) -> tuple[str, str, str, str]:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))
    table_path = os.path.join(output_dir, "integration_comparison_table.txt")
    csv_path = os.path.join(output_dir, "integration_comparison_results.csv")
    analysis_path = os.path.join(output_dir, "integration_comparison_analysis.txt")
    graph_path = os.path.join(output_dir, "integration_comparison_graph.png")

    columns = ["h", "method", "n", "result", "exact", "absolute_error"]
    with open(table_path, "w", encoding="utf-8") as file:
        file.write("Integration Method Comparison\n\n")
        file.write(
            f"{'h':<14}{'Method':<30}{'Points':<10}"
            f"{'Approximation':<20}{'Exact':<20}{'Absolute Error':<20}\n"
        )
        for row in rows:
            exact = "N/A" if row["exact"] is None else f"{row['exact']:.12g}"
            error = (
                "N/A"
                if row["absolute_error"] is None
                else f"{row['absolute_error']:.12g}"
            )
            file.write(
                f"{row['h']:<14.6g}{row['method']:<30}{row['n']:<10}"
                f"{row['result']:<20.12g}{exact:<20}{error:<20}\n"
            )

    with open(csv_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[float, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(float(row["h"]), []).append(row)
    with open(analysis_path, "w", encoding="utf-8") as file:
        file.write("Integration Method Comparison Analysis\n")
        for h, results in sorted(grouped.items(), reverse=True):
            file.write(f"\nh = {h:.12g}\n")
            errors = [row for row in results if row["absolute_error"] is not None]
            if errors:
                best = min(errors, key=lambda row: float(row["absolute_error"]))
                file.write(
                    f"Most accurate: {best['method']} "
                    f"(absolute error {best['absolute_error']:.12g})\n"
                )
            for row in results:
                file.write(f"{row['method']}: {row['result']:.12g}\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    methods = (
        ("CompositeTrapezoidal", "Composite Trapezoidal"),
        ("CompositeSimpsonOneThird", "Composite Simpson's 1/3"),
        ("CompositeSimpsonThreeEighths", "Composite Simpson's 3/8"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
    for axis, (method_name, title) in zip(axes, methods):
        method_rows = [
            row
            for row in rows
            if row["method"] == method_name and row["absolute_error"] is not None
        ]
        if not method_rows:
            continue

        h_values = [float(row["h"]) for row in method_rows]
        errors = [max(float(row["absolute_error"]), 1e-16) for row in method_rows]
        axis.loglog(
            h_values,
            errors,
            marker="o",
            linestyle="-",
            label=title,
        )
        axis.set_title(title)
        axis.set_xlabel("h")
        axis.set_ylabel("Absolute error")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)

    axes[0].set_ylabel("Absolute error")
    figure.suptitle("Integration Comparison Error vs Step Size", fontsize=15)
    figure.tight_layout()
    figure.savefig(graph_path, dpi=180)
    plt.close(figure)
    return table_path, csv_path, analysis_path, graph_path


def calculate_integration_comparison(
    data_source, exact: float | None, h_values: list[float]
) -> list[dict[str, object]]:
    method_types = (
        CompositeTrapezoidal,
        CompositeSimpsonOneThird,
        CompositeSimpsonThreeEighths,
    )
    rows: list[dict[str, object]] = []
    for h in h_values:
        if isinstance(data_source, FunctionDataSource):
            intervals = round(abs(data_source.b - data_source.a) / h)
            current_source = FunctionDataSource(
                data_source.func, data_source.a, data_source.b, intervals + 1
            )
            actual_h = abs(current_source.x_values[1] - current_source.x_values[0])
        else:
            current_source = data_source
            actual_h = h

        for method_type in method_types:
            result = method_type(current_source).result_row(exact=exact)
            rows.append(
                {
                    "h": actual_h,
                    "method": result.method,
                    "n": result.n,
                    "result": result.result,
                    "exact": result.exact,
                    "absolute_error": result.absolute_error,
                }
            )
    return rows


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


def run_integration() -> None:
    print_line("NUMERICAL INTEGRATION")
    input_path = os.path.join(PROJECT_DIR, "input", "input_integration.txt")
    try:
        data_source, exact, h_values = read_integration_input(input_path)
    except (FileNotFoundError, ValueError) as error:
        print(f"Input error: {error}")
        return

    rows = calculate_integration_comparison(data_source, exact, h_values)
    print(f"Step sizes: {h_values}")
    print("Compared methods: trapezoidal, Simpson's 1/3, Simpson's 3/8")
    output_paths = save_integration_comparison(rows, data_source)
    print("\nComparison results saved:")
    for path in output_paths:
        print(os.path.relpath(path, PROJECT_DIR))


def run_eigenvalues() -> None:
    print_line("EIGENVALUES AND EIGENVECTORS")
    print("Choose a solver:")
    print("1. Power Method         (input/input_eigen_power.txt)")
    print("2. Inverse Power Method (input/input_eigen_inverse.txt)")
    print("3. Shifted Inverse Power Method (input/input_eigen_shifted.txt)")
    print("4. Jacobi Method        (input/input_eigen_jacobi.txt)")

    choice = input("Enter your choice: ").strip()
    method_files = {
        "1": "input_eigen_power.txt",
        "2": "input_eigen_inverse.txt",
        "3": "input_eigen_shifted.txt",
        "4": "input_eigen_jacobi.txt",
    }
    file_name = method_files.get(choice)
    if file_name is None:
        print("Invalid choice.")
        return

    input_path = os.path.join(PROJECT_DIR, "input", file_name)
    try:
        config = read_eigen_input(input_path)
    except (FileNotFoundError, ValueError) as error:
        print(f"Input error: {error}")
        return

    method_name = config["method"]
    matrix = config["matrix"]
    tolerance = config.get("tolerance", 1e-10)
    max_iterations = config.get("max_iterations", 1000)

    if method_name == "power":
        solver = PowerMethod(matrix, tolerance=tolerance, max_iterations=max_iterations, initial_vector=config.get("initial_vector"))
        result = solver.solve()
    elif method_name == "inverse":
        solver = InversePowerMethod(matrix, tolerance=tolerance, max_iterations=max_iterations, initial_vector=config.get("initial_vector"))
        result = solver.solve()
    elif method_name == "shifted":
        solver = ShiftedInversePowerMethod(
            matrix,
            mu=config["mu"],
            tolerance=tolerance,
            max_iterations=max_iterations,
            initial_vector=config.get("initial_vector"),
        )
        result = solver.solve()
    elif method_name == "jacobi":
        solver = JacobiMethod(matrix, tolerance=tolerance, max_iterations=max_iterations)
        result = solver.solve()
    else:
        raise ValueError(f"Unsupported eigen method: {method_name}")

    path = os.path.join(PROJECT_DIR, "output", f"{method_name}_results.txt")
    save_eigen_results(result, path)
    print("\nEigenvalue computation complete.")
    if "eigenvalue" in result:
        print(f"Eigenvalue: {result['eigenvalue']:.12g}")
        print(f"Residual norm: {result['residual_norm']:.12g}")
    else:
        print(f"Eigenvalues: {[float(v) for v in result['eigenvalues']]}")
        print(f"Residual norms: {[float(v) for v in result['residual_norms']]}")
    print(f"Results saved to {os.path.relpath(path, PROJECT_DIR)}")


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
        print("4. Numerical Integration")
        print("5. Eigenvalues and Eigenvectors")
        print("0. Exit")

        choice = input("Enter your choice: ").strip()

        if choice == "1":
            run_differentiation()
        elif choice == "2":
            run_richardson()
        elif choice == "3":
            run_interpolation()
        elif choice == "4":
            run_integration()
        elif choice == "5":
            run_eigenvalues()
        elif choice == "0":
            print("Goodbye.")
            break
        else:
            print("Invalid choice.")


if __name__ == "__main__":
    main()
