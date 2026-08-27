import math
import os
import sys

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(PROJECT_DIR, "src")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from analyzer import DifferentiationAnalyzer, RichardsonAnalyzer
from utils import ensure_output_dir, read_step_size_input


TEST_FUNCTIONS = {
    "1": ("exp", math.exp, math.exp),
    "2": ("sin", math.sin, math.cos),
    "3": ("cos", math.cos, lambda x: -math.sin(x)),
    "4": (
        "poly",
        lambda x: x ** 3 - 2 * x + 1,
        lambda x: 3 * x ** 2 - 2,
    ),
}


def print_line(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def choose_functions():
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


def read_input_file(filename: str):
    path = os.path.join(PROJECT_DIR, "input", filename)

    try:
        return read_step_size_input(path)
    except (FileNotFoundError, ValueError) as error:
        print(f"Error: {error}")
        return None


def save_differentiation_results(analyzer: DifferentiationAnalyzer) -> None:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))

    analyzer.write_table(
        os.path.join(output_dir, "differentiation_table.txt")
    )
    analyzer.write_csv(
        os.path.join(output_dir, "differentiation_results.csv")
    )
    analyzer.write_analysis(
        os.path.join(output_dir, "differentiation_analysis.txt")
    )
    analyzer.plot(
        os.path.join(output_dir, "differentiation_loglog_plot.png")
    )


def save_richardson_results(analyzer: RichardsonAnalyzer) -> None:
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))

    analyzer.write_table(
        os.path.join(output_dir, "richardson_table.txt")
    )
    analyzer.write_csv(
        os.path.join(output_dir, "richardson_results.csv")
    )
    analyzer.write_analysis(
        os.path.join(output_dir, "richardson_analysis.txt")
    )
    analyzer.plot(
        os.path.join(output_dir, "richardson_loglog_plot.png")
    )


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

    for name, function, exact_derivative in functions:
        analyzer.add_function(name, function, exact_derivative)

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

    for name, function, exact_derivative in functions:
        analyzer.add_function(name, function, exact_derivative)

    analyzer.run()

    print("\n" + analyzer.format_table())
    print(analyzer.format_analysis())

    save_richardson_results(analyzer)
    print("\nResults saved in the output folder.")


def main() -> None:
    while True:
        print_line("NUMERICAL COMPUTING PROJECT")
        print("1. Numerical Differentiation")
        print("2. Richardson Extrapolation")
        print("0. Exit")

        choice = input("Enter your choice: ").strip()

        if choice == "1":
            run_differentiation()
        elif choice == "2":
            run_richardson()
        elif choice == "0":
            print("Goodbye.")
            break
        else:
            print("Invalid choice.")


if __name__ == "__main__":
    main()
