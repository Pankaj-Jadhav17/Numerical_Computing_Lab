import csv
import math
from dataclasses import asdict, dataclass
from typing import Callable, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from differentiation import (
    BackwardDifference,
    CentralDifference,
    ForwardDifference,
)


@dataclass
class DiffResultRow:
    function: str
    h: float
    forward: float
    backward: float
    central: float
    exact: float
    err_forward: float
    err_backward: float
    err_central: float


class DifferentiationAnalyzer:
    def __init__(self, x0: float, h_values: List[float]):
        self.x0 = x0
        self.h_values = h_values
        self.test_functions: List[Dict] = []
        self.results: List[DiffResultRow] = []

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
                forward = ForwardDifference(test["function"], h).derivative(self.x0)
                backward = BackwardDifference(test["function"], h).derivative(self.x0)
                central = CentralDifference(test["function"], h).derivative(self.x0)

                self.results.append(
                    DiffResultRow(
                        function=test["name"],
                        h=h,
                        forward=forward,
                        backward=backward,
                        central=central,
                        exact=exact,
                        err_forward=abs(exact - forward),
                        err_backward=abs(exact - backward),
                        err_central=abs(exact - central),
                    )
                )

    def grouped_results(self) -> Dict[str, List[DiffResultRow]]:
        groups: Dict[str, List[DiffResultRow]] = {}

        for result in self.results:
            groups.setdefault(result.function, []).append(result)

        return groups

    def format_table(self) -> str:
        if not self.results:
            return "No results available."

        header = (
            f"{'Function':<10}{'h':<12}{'Forward':<16}"
            f"{'Backward':<16}{'Central':<16}{'Exact':<16}"
            f"{'Err Forward':<16}{'Err Backward':<16}{'Err Central':<16}\n"
        )

        lines = [header]

        for result in self.results:
            lines.append(
                f"{result.function:<10}"
                f"{result.h:<12.0e}"
                f"{result.forward:<16.6e}"
                f"{result.backward:<16.6e}"
                f"{result.central:<16.6e}"
                f"{result.exact:<16.6e}"
                f"{result.err_forward:<16.6e}"
                f"{result.err_backward:<16.6e}"
                f"{result.err_central:<16.6e}\n"
            )

        return "".join(lines)

    def write_table(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as file:
            file.write("Numerical Differentiation Results\n")
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

        columns = min(3, len(groups))
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
                self._plot_errors([item.err_forward for item in results]),
                marker="o",
                label="Forward",
            )
            axis.loglog(
                h_values,
                self._plot_errors([item.err_backward for item in results]),
                marker="s",
                label="Backward",
            )
            axis.loglog(
                h_values,
                self._plot_errors([item.err_central for item in results]),
                marker="^",
                label="Central",
            )

            axis.set_title(f"f(x) = {name}")
            axis.set_xlabel("h")
            axis.set_ylabel("Absolute Error")
            axis.invert_xaxis()
            axis.grid(True, which="both")
            axis.legend()

        for axis in axes[len(groups):]:
            axis.axis("off")

        figure.suptitle("Numerical Differentiation Error")
        figure.tight_layout()
        figure.savefig(path, dpi=200)
        plt.close(figure)

    @staticmethod
    def _plot_errors(errors: List[float]) -> List[float]:
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

        lines = ["\n===== ANALYSIS =====\n"]

        for name, results in self.grouped_results().items():
            middle = results[len(results) // 2]

            lines.append(f"\nFunction: {name}\n")
            lines.append(
                f"At h = {middle.h:.0e}: "
                f"Forward Error = {middle.err_forward:.3e}, "
                f"Backward Error = {middle.err_backward:.3e}, "
                f"Central Error = {middle.err_central:.3e}\n"
            )

            errors = {
                "Forward": middle.err_forward,
                "Backward": middle.err_backward,
                "Central": middle.err_central,
            }

            best_method = min(errors, key=errors.get)
            lines.append(f"Most accurate method: {best_method}\n")

            if len(results) >= 2:
                forward_order = self.observed_order(
                    results[0].h,
                    results[0].err_forward,
                    results[1].h,
                    results[1].err_forward,
                )
                backward_order = self.observed_order(
                    results[0].h,
                    results[0].err_backward,
                    results[1].h,
                    results[1].err_backward,
                )
                central_order = self.observed_order(
                    results[0].h,
                    results[0].err_central,
                    results[1].h,
                    results[1].err_central,
                )

                lines.append(
                    f"Observed order: Forward = {forward_order:.3f}, "
                    f"Backward = {backward_order:.3f}, "
                    f"Central = {central_order:.3f}\n"
                )

            minimum = min(results, key=lambda item: item.err_central)
            lines.append(
                f"Minimum Central Difference error: "
                f"{minimum.err_central:.3e} at h = {minimum.h:.0e}\n"
            )

        return "".join(lines)

    def write_analysis(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as file:
            file.write("Numerical Differentiation Analysis\n")
            file.write(f"Evaluation point: {self.x0}\n")
            file.write(self.format_analysis())
