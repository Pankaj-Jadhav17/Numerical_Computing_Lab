#!/usr/bin/env python3
"""Quickstart script: run analyzers non-interactively and save outputs.

Usage: python quickstart.py
"""
import os
import sys
import math

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(PROJECT_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from analyzer import DifferentiationAnalyzer, InterpolationAnalyzer, RichardsonAnalyzer
from utils.io_utils import read_step_size_input, read_interpolation_input, ensure_output_dir


TEST_FUNCTIONS = [
    ("exp", math.exp, math.exp),
    ("sin", math.sin, math.cos),
    ("cos", math.cos, lambda x: -math.sin(x)),
    ("poly", (lambda x: x**3 - 2 * x + 1), (lambda x: 3 * x**2 - 2)),
]


def _save_analysis_outputs(
    analyzer,
    output_dir,
    table_name,
    csv_name,
    analysis_name,
    plot_name,
):
    analyzer.write_table(os.path.join(output_dir, table_name))
    analyzer.write_csv(os.path.join(output_dir, csv_name))
    analyzer.write_analysis(os.path.join(output_dir, analysis_name))
    if hasattr(analyzer, "plot"):
        analyzer.plot(os.path.join(output_dir, plot_name))


def run_all():
    output_dir = ensure_output_dir(os.path.join(PROJECT_DIR, "output"))

    # Differentiation
    try:
        x0, h_values = read_step_size_input(os.path.join(PROJECT_DIR, "input", "input_differentiation.txt"))
        analyzer = DifferentiationAnalyzer(x0, h_values)
        for name, func, exact in TEST_FUNCTIONS:
            analyzer.add_function(name, func, exact)
        analyzer.run()
        _save_analysis_outputs(
            analyzer,
            output_dir,
            "differentiation_table.txt",
            "differentiation_results.csv",
            "differentiation_analysis.txt",
            "differentiation_loglog_plot.png",
        )
        print("Differentiation analysis complete.")
    except Exception as exc:
        print(f"Skipping differentiation: {exc}")

    # Richardson
    try:
        x0, h_values = read_step_size_input(os.path.join(PROJECT_DIR, "input", "input_richardson.txt"))
        analyzer = RichardsonAnalyzer(x0, h_values)
        for name, func, exact in TEST_FUNCTIONS:
            analyzer.add_function(name, func, exact)
        analyzer.run()
        _save_analysis_outputs(
            analyzer,
            output_dir,
            "richardson_table.txt",
            "richardson_results.csv",
            "richardson_analysis.txt",
            "richardson_loglog_plot.png",
        )
        print("Richardson analysis complete.")
    except Exception as exc:
        print(f"Skipping Richardson: {exc}")

    # Interpolation
    try:
        x_data, y_data, query_points = read_interpolation_input(os.path.join(PROJECT_DIR, "input", "input_interpolation.txt"))
        interp = InterpolationAnalyzer(query_points)
        for name, func, _ in TEST_FUNCTIONS:
            interp.add_function(name, func, x_data, y_data if y_data else None)
        interp.run()
        _save_analysis_outputs(
            interp,
            output_dir,
            "interpolation_table.txt",
            "interpolation_results.csv",
            "interpolation_analysis.txt",
            "interpolation_plot.png",
        )
        print("Interpolation analysis complete.")
    except Exception as exc:
        print(f"Skipping interpolation: {exc}")


if __name__ == "__main__":
    run_all()
