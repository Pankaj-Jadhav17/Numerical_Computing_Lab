"""
InterpolationAnalyzer — orchestrates Newton interpolation runs and
generates formatted tables, CSV exports, analysis text, and plots.

Follows the same pattern as DifferentiationAnalyzer and RichardsonAnalyzer.
"""

import csv
import math
from collections.abc import Callable
from dataclasses import asdict, dataclass
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from interpolation import NewtonInterpolation


@dataclass
class InterpolationResultRow:
    """One row of interpolation results for a single query point."""
    function_name: str
    n_points: int
    x_query: float
    interpolated: float
    exact: float
    abs_error: float
    error_estimate: float


class InterpolationAnalyzer:
    """Analyzer for Newton Interpolation with Divided Differences.

    Workflow:
        1. Create with query points.
        2. add_function(...) to register test functions.
        3. run() to compute results for each function on small & large datasets.
        4. Use format_table(), format_analysis(), plot(), write_*() for output.
    """

    def __init__(self, query_points: List[float]) -> None:
        self.query_points = query_points
        self.test_functions: list[dict] = []
        self.results: List[InterpolationResultRow] = []
        self._interpolators: dict[str, dict[int, NewtonInterpolation]] = {}

    def add_function(
        self,
        name: str,
        function: Callable[[float], float],
        x_data: List[float],
        y_data: Optional[List[float]] = None,
    ) -> None:
        """Register a function to interpolate.

        If y_data is None, it is computed from function(x_data[i]).
        """
        if y_data is None:
            y_data = [function(x) for x in x_data]

        self.test_functions.append({
            "name": name,
            "function": function,
            "x_data": list(x_data),
            "y_data": list(y_data),
        })

    def run(self) -> None:
        """Run Newton interpolation for every registered function and
        query point.  For each function the interpolation is tested at
        multiple dataset sizes: full, half, and quarter (when possible)."""
        self.results.clear()
        self._interpolators.clear()

        for test in self.test_functions:
            name = test["name"]
            func = test["function"]
            x_full = test["x_data"]
            y_full = test["y_data"]
            n_full = len(x_full)

            # Dataset sizes to test: quarter, half, full
            sizes = sorted({
                max(2, n_full // 4),
                max(2, n_full // 2),
                n_full,
            })

            self._interpolators[name] = {}

            for size in sizes:
                x_sub = x_full[:size]
                y_sub = y_full[:size]

                interp = NewtonInterpolation(x_sub, y_sub)
                self._interpolators[name][size] = interp

                for xq in self.query_points:
                    value, error_est = interp.evaluate_with_error(xq)
                    exact = func(xq)

                    self.results.append(InterpolationResultRow(
                        function_name=name,
                        n_points=size,
                        x_query=xq,
                        interpolated=value,
                        exact=exact,
                        abs_error=abs(exact - value),
                        error_estimate=error_est,
                    ))

    # ── grouped access ─────────────────────────────────────────────────
    def grouped_results(self) -> dict[str, List[InterpolationResultRow]]:
        groups: dict[str, List[InterpolationResultRow]] = {}
        for r in self.results:
            groups.setdefault(r.function_name, []).append(r)
        return groups

    # ── formatted output ───────────────────────────────────────────────
    def format_table(self) -> str:
        if not self.results:
            return "No results available."

        header = (
            f"{'Function':<10}{'n':<6}{'x_query':<12}"
            f"{'Interpolated':<18}{'Exact':<18}"
            f"{'Abs Error':<16}{'Error Est':<16}\n"
        )
        lines = [header, "-" * len(header.rstrip()) + "\n"]

        for r in self.results:
            lines.append(
                f"{r.function_name:<10}"
                f"{r.n_points:<6}"
                f"{r.x_query:<12.4f}"
                f"{r.interpolated:<18.8e}"
                f"{r.exact:<18.8e}"
                f"{r.abs_error:<16.8e}"
                f"{r.error_estimate:<16.8e}\n"
            )
        return "".join(lines)

    def format_dd_tables(self) -> str:
        """Return the divided-difference tables for all registered functions."""
        if not self._interpolators:
            return "Run the analyzer first."

        lines: List[str] = []
        for name, size_map in self._interpolators.items():
            for size, interp in size_map.items():
                lines.append(f"\n{'='*60}")
                lines.append(f"  Divided Difference Table: {name} (n={size})")
                lines.append(f"{'='*60}")
                lines.append(interp.format_table())
                lines.append("")
        return "\n".join(lines)

    def write_table(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write("Newton Interpolation Results\n")
            f.write(f"Query points: {self.query_points}\n\n")
            f.write(self.format_table())
            f.write("\n\n")
            f.write(self.format_dd_tables())

    def write_csv(self, path: str) -> None:
        if not self.results:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            fieldnames = list(asdict(self.results[0]).keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in self.results:
                writer.writerow(asdict(r))

    # ── plotting ───────────────────────────────────────────────────────
    def plot(self, path: str) -> None:
        """Generate a multi-panel plot:
           - Top row: interpolation curves vs exact function
           - Bottom row: error vs number of data points
        """
        groups = self.grouped_results()
        if not groups:
            return

        n_funcs = len(groups)
        fig, axes = plt.subplots(
            2, n_funcs,
            figsize=(7 * n_funcs, 10),
            squeeze=False,
        )

        for col, (name, results) in enumerate(groups.items()):
            ax_top = axes[0][col]
            ax_bot = axes[1][col]

            # ── top panel: curves ──────────────────────────────────────
            interp_map = self._interpolators.get(name, {})
            if interp_map:
                # Find x range
                all_x = []
                for interp in interp_map.values():
                    all_x.extend(interp.x_data)
                x_min, x_max = min(all_x), max(all_x)
                margin = (x_max - x_min) * 0.05
                x_plot = [
                    x_min - margin + i * (x_max - x_min + 2 * margin) / 200
                    for i in range(201)
                ]

                # Exact curve
                func = None
                for t in self.test_functions:
                    if t["name"] == name:
                        func = t["function"]
                        break

                if func:
                    y_exact = [func(x) for x in x_plot]
                    ax_top.plot(x_plot, y_exact, "k-", lw=2, label="Exact")

                # Interpolated curves for each dataset size
                colors = plt.cm.viridis([0.2, 0.5, 0.8])
                for idx, (size, interp) in enumerate(
                    sorted(interp_map.items())
                ):
                    y_interp = interp.interpolate_many(x_plot)
                    ax_top.plot(
                        x_plot, y_interp,
                        "--",
                        color=colors[idx % len(colors)],
                        lw=1.5,
                        label=f"n={size}",
                    )
                    # Plot data points
                    ax_top.plot(
                        interp.x_data, interp.y_data,
                        "o",
                        color=colors[idx % len(colors)],
                        markersize=5,
                    )

            ax_top.set_title(f"f(x) = {name}")
            ax_top.set_xlabel("x")
            ax_top.set_ylabel("f(x)")
            ax_top.legend(fontsize=8)
            ax_top.grid(True, alpha=0.3)

            # ── bottom panel: error vs n ───────────────────────────────
            # Group results by query point
            by_query: dict[float, list[InterpolationResultRow]] = {}
            for r in results:
                by_query.setdefault(r.x_query, []).append(r)

            for xq, rows in by_query.items():
                rows_sorted = sorted(rows, key=lambda r: r.n_points)
                ns = [r.n_points for r in rows_sorted]
                errs = [max(r.abs_error, 1e-16) for r in rows_sorted]
                ax_bot.semilogy(ns, errs, "o-", label=f"x={xq:.2f}")

            ax_bot.set_title(f"Error vs Dataset Size: {name}")
            ax_bot.set_xlabel("Number of data points (n)")
            ax_bot.set_ylabel("Absolute Error")
            ax_bot.legend(fontsize=8)
            ax_bot.grid(True, which="both", alpha=0.3)

        fig.suptitle(
            "Newton Interpolation with Divided Differences",
            fontsize=14, fontweight="bold",
        )
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)

    # ── analysis ───────────────────────────────────────────────────────
    def format_analysis(self) -> str:
        if not self.results:
            return "No analysis available."

        lines = ["\n===== NEWTON INTERPOLATION ANALYSIS =====\n"]

        for name, results in self.grouped_results().items():
            lines.append(f"\nFunction: {name}")
            lines.append(f"  Total evaluations: {len(results)}")

            # Group by dataset size
            by_size: dict[int, list[InterpolationResultRow]] = {}
            for r in results:
                by_size.setdefault(r.n_points, []).append(r)

            for size in sorted(by_size.keys()):
                rows = by_size[size]
                max_err = max(r.abs_error for r in rows)
                avg_err = sum(r.abs_error for r in rows) / len(rows)
                min_err = min(r.abs_error for r in rows)

                lines.append(f"\n  Dataset size n = {size}:")
                lines.append(f"    Max absolute error:  {max_err:.6e}")
                lines.append(f"    Avg absolute error:  {avg_err:.6e}")
                lines.append(f"    Min absolute error:  {min_err:.6e}")

            # Convergence: compare smallest vs largest dataset
            sizes = sorted(by_size.keys())
            if len(sizes) >= 2:
                small_avg = sum(
                    r.abs_error for r in by_size[sizes[0]]
                ) / len(by_size[sizes[0]])
                large_avg = sum(
                    r.abs_error for r in by_size[sizes[-1]]
                ) / len(by_size[sizes[-1]])

                if small_avg > 0 and large_avg > 0:
                    improvement = small_avg / large_avg
                    lines.append(
                        f"\n  Accuracy improvement (n={sizes[0]} → n={sizes[-1]}): "
                        f"{improvement:.1f}x"
                    )

            lines.append("")

        lines.append("\n===== NOTES =====")
        lines.append(
            "• Newton interpolation with n points produces a polynomial of degree n-1."
        )
        lines.append(
            "• Adding more data points generally improves accuracy inside the data range."
        )
        lines.append(
            "• Extrapolation outside the data range can produce large errors."
        )
        lines.append(
            "• The error estimate uses the last divided-difference coefficient "
            "as an approximation."
        )
        lines.append("")

        return "\n".join(lines)

    def write_analysis(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write("Newton Interpolation Analysis\n")
            f.write(f"Query points: {self.query_points}\n")
            f.write(self.format_analysis())
