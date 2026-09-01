import csv
import math
from collections.abc import Callable
from dataclasses import asdict, dataclass

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

from differentiation import CentralDifference, RichardsonExtrapolation


@dataclass
class RichardsonResultRow:
    function: str
    h: float
    central: float
    richardson: float
    exact: float
    err_central: float
    err_richardson: float


class RichardsonAnalyzer:
    def __init__(self, x0: float, h_values: list[float]):
        self.x0 = x0
        self.h_values = h_values
        self.test_functions: list[dict[str, object]] = []
        self.results: list[RichardsonResultRow] = []

    def add_function(
        self,
        name: str,
        function: Callable[[float], float],
        exact_derivative: Callable[[float], float],
    ) -> None:
        self.test_functions.append(
            {
                "name": name,
                "function": function,
                "exact": exact_derivative,
            }
        )

    def run(self) -> None:
        self.results.clear()

        for test in self.test_functions:
            exact = test["exact"](self.x0)

            for h in self.h_values:
                central = CentralDifference(test["function"], h).derivative(self.x0)
                richardson = RichardsonExtrapolation(
                    test["function"], h
                ).derivative(self.x0)

                self.results.append(
                    RichardsonResultRow(
                        function=test["name"],
                        h=h,
                        central=central,
                        richardson=richardson,
                        exact=exact,
                        err_central=abs(exact - central),
                        err_richardson=abs(exact - richardson),
                    )
                )

    def grouped_results(self) -> dict[str, list[RichardsonResultRow]]:
        groups: dict[str, list[RichardsonResultRow]] = {}

        for result in self.results:
            groups.setdefault(result.function, []).append(result)

        return groups

    def format_table(self) -> str:
        if not self.results:
            return "No results available."

        header = (
            f"{'Function':<10}{'h':<12}{'Central':<18}"
            f"{'Richardson':<18}{'Exact':<18}"
            f"{'Err Central':<18}{'Err Richardson':<18}\n"
        )

        lines = [header]

        for result in self.results:
            lines.append(
                f"{result.function:<10}"
                f"{result.h:<12.0e}"
                f"{result.central:<18.8e}"
                f"{result.richardson:<18.8e}"
                f"{result.exact:<18.8e}"
                f"{result.err_central:<18.8e}"
                f"{result.err_richardson:<18.8e}\n"
            )

        return "".join(lines)

    def write_table(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as file:
            file.write("Richardson Extrapolation Results\n")
            file.write(f"Evaluation point: {self.x0}\n\n")
            file.write(self.format_table())

    def write_csv(self, path: str) -> None:
        if not self.results:
            return

        with open(path, "w", newline="", encoding="utf-8") as file:
            fieldnames = list(asdict(self.results[0]).keys())
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()

            for result in self.results:
                writer.writerow(asdict(result))

    def plot(self, path: str) -> None:
        groups = self.grouped_results()

        if not groups:
            return

        columns = min(2, len(groups))
        rows = math.ceil(len(groups) / columns)

        figure, axes = plt.subplots(
            rows,
            columns,
            figsize=(6 * columns, 4.5 * rows),
            squeeze=False,
        )

        axes = axes.flatten()

        for axis, (name, results) in zip(axes, groups.items()):
            h_values = [item.h for item in results]

            axis.loglog(
                h_values,
                self._plot_errors([item.err_central for item in results]),
                marker="o",
                label="Central O(h^2)",
            )
            axis.loglog(
                h_values,
                self._plot_errors([item.err_richardson for item in results]),
                marker="^",
                label="Richardson O(h^4)",
            )

            axis.set_title(f"f(x) = {name}")
            axis.set_xlabel("h")
            axis.set_ylabel("Absolute Error")
            axis.invert_xaxis()
            axis.grid(True, which="both")
            axis.legend()

        for axis in axes[len(groups):]:
            axis.axis("off")

        figure.suptitle("Central Difference vs Richardson Extrapolation")
        figure.tight_layout()
        figure.savefig(path, dpi=200)
        plt.close(figure)

    @staticmethod
    def _plot_errors(errors: list[float]) -> list[float]:
        smallest = 1e-16
        return [max(error, smallest) for error in errors]

    @staticmethod
    def observed_order(
        h1: float,
        error1: float,
        h2: float,
        error2: float,
    ) -> float:
        if error1 <= 0 or error2 <= 0 or h1 == h2:
            return 0.0

        return math.log(error1 / error2) / math.log(h1 / h2)

    def format_analysis(self) -> str:
        if not self.results:
            return "No analysis available."

        lines = ["\n===== RICHARDSON ANALYSIS =====\n"]

        for name, results in self.grouped_results().items():
            lines.append(f"\nFunction: {name}\n")

            if len(results) >= 2:
                central_order = self.observed_order(
                    results[0].h,
                    results[0].err_central,
                    results[1].h,
                    results[1].err_central,
                )
                richardson_order = self.observed_order(
                    results[0].h,
                    results[0].err_richardson,
                    results[1].h,
                    results[1].err_richardson,
                )

                lines.append(
                    f"Observed order: Central = {central_order:.3f}, "
                    f"Richardson = {richardson_order:.3f}\n"
                )

            better_count = sum(
                item.err_richardson < item.err_central for item in results
            )

            lines.append(
                f"Richardson is more accurate for "
                f"{better_count} out of {len(results)} tested h values.\n"
            )

            central_minimum = min(results, key=lambda item: item.err_central)
            richardson_minimum = min(
                results, key=lambda item: item.err_richardson
            )

            lines.append(
                f"Minimum Central error: {central_minimum.err_central:.3e} "
                f"at h = {central_minimum.h:.0e}\n"
            )
            lines.append(
                f"Minimum Richardson error: {richardson_minimum.err_richardson:.3e} "
                f"at h = {richardson_minimum.h:.0e}\n"
            )

        lines.append("\n===== ASSIGNMENT QUESTIONS =====\n")
        lines.append(
            "Q1. Does Richardson always give a smaller error than Central Difference?\n"
        )
        lines.append(
            "Richardson usually gives a smaller truncation error, but it may not "
            "always be better when h becomes very small because round-off error can dominate.\n\n"
        )
        lines.append(
            "Q2. What happens when h is reduced by 10?\n"
        )
        lines.append(
            "For Central Difference, the theoretical error decreases by about 100 times. "
            "For Richardson, the theoretical error decreases by about 10000 times.\n\n"
        )
        lines.append(
            "Q3. Does Richardson continue to improve for very small h?\n"
        )
        lines.append(
            "No. Floating-point round-off error eventually becomes important and the total error can increase.\n\n"
        )
        lines.append(
            "Q4. Does the log-log slope support the theory?\n"
        )
        lines.append(
            "The slope should be close to 2 for Central Difference and close to 4 for Richardson "
            "before round-off error becomes dominant.\n\n"
        )
        lines.append(
            "Q5. Does Richardson achieve O(h^4)?\n"
        )
        lines.append(
            "Yes, when truncation error is the main source of error. Richardson removes the leading O(h^2) term.\n\n"
        )
        lines.append(
            "Q6. What is the role of floating-point round-off error?\n"
        )
        lines.append(
            "When h becomes very small, subtracting nearly equal floating-point values loses accuracy. "
            "Dividing by a small h increases the effect of this error.\n"
        )

        return "".join(lines)

    def write_analysis(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as file:
            file.write("Richardson Extrapolation Analysis\n")
            file.write(f"Evaluation point: {self.x0}\n")
            file.write(self.format_analysis())
