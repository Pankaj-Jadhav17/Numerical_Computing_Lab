import matplotlib.pyplot as plt
import pytest

from integration import (
    ArrayDataSource,
    CompositeSimpsonOneThird,
    CompositeSimpsonThreeEighths,
    CompositeTrapezoidal,
    FunctionDataSource,
    IntegrationMethod,
)
import main
from main import (
    calculate_integration_comparison,
    save_integration_comparison,
    save_integration_result,
)
from utils import read_integration_input


@pytest.mark.parametrize(
    ("method_type", "a", "b", "n"),
    [
        (CompositeSimpsonOneThird, 0.0, 4.0, 5),
        (CompositeSimpsonThreeEighths, 0.0, 6.0, 7),
    ],
)
def test_simpson_methods_integrate_cubic_exactly(method_type, a, b, n):
    source = FunctionDataSource(lambda x: x**3 - 2 * x + 1, a, b, n)
    result = method_type(source).integrate()
    exact = (b**4 / 4 - b**2 + b) - (a**4 / 4 - a**2 + a)

    assert result == pytest.approx(exact)


@pytest.mark.parametrize(
    ("method_type", "n"),
    [(CompositeSimpsonOneThird, 4), (CompositeSimpsonThreeEighths, 5)],
)
def test_simpson_methods_reject_invalid_subinterval_count(method_type, n):
    source = FunctionDataSource(lambda x: x**2, 0.0, 1.0, n)

    with pytest.raises(ValueError, match="subintervals"):
        method_type(source).integrate()


@pytest.mark.parametrize(
    ("method_type", "x_values"),
    [
        (CompositeSimpsonOneThird, [0.0, 1.0, 2.1]),
        (CompositeSimpsonThreeEighths, [0.0, 1.0, 2.0, 3.1]),
    ],
)
def test_simpson_methods_reject_nonuniform_points(method_type, x_values):
    source = ArrayDataSource(x_values, [x**2 for x in x_values])

    with pytest.raises(ValueError, match="equally spaced"):
        method_type(source).integrate()


@pytest.mark.parametrize(
    ("method_type", "n"),
    [(CompositeSimpsonOneThird, 5), (CompositeSimpsonThreeEighths, 4)],
)
def test_simpson_methods_handle_reversed_interval_and_result_row(method_type, n):
    source = FunctionDataSource(lambda x: x**3, 3.0, 0.0, n)
    method = method_type(source)

    assert isinstance(method, IntegrationMethod)
    assert method.integrate() == pytest.approx(-81 / 4)
    row = method.result_row(exact=-81 / 4)
    assert row.result == pytest.approx(-81 / 4)
    assert row.absolute_error == pytest.approx(0.0)


def test_read_integration_input_for_function(tmp_path):
    input_file = tmp_path / "integration.txt"
    input_file.write_text(
        "source=function\nfunction=poly\na=0\nb=4\n"
        "h=0.6666666666667, 0.3333333333333\nexact=52\n",
        encoding="utf-8",
    )

    data_source, exact, h_values = read_integration_input(str(input_file))

    assert data_source.x_values == pytest.approx([index * 2 / 3 for index in range(7)])
    assert data_source.y_values[0] == pytest.approx(1)
    assert data_source.y_values[-1] == pytest.approx(57)
    assert exact == 52
    assert h_values == pytest.approx([2 / 3, 1 / 3])


def test_read_integration_input_for_arrays(tmp_path):
    input_file = tmp_path / "integration.txt"
    x_values = list(range(7))
    y_values = [x**2 for x in x_values]
    input_file.write_text(
        "source=array\nx=" + ", ".join(map(str, x_values))
        + "\ny=" + ", ".join(map(str, y_values)),
        encoding="utf-8",
    )

    data_source, exact, h_values = read_integration_input(str(input_file))

    assert data_source.x_values == x_values
    assert data_source.y_values == y_values
    assert exact is None
    assert h_values == [1.0]


def test_read_integration_input_rejects_step_size_without_shared_panel_count(tmp_path):
    input_file = tmp_path / "integration.txt"
    input_file.write_text(
        "source=function\nfunction=sin\na=0\nb=1\nh=0.2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="divisible by 6"):
        read_integration_input(str(input_file))


def test_calculate_integration_comparison_uses_all_rules_for_each_h():
    source = FunctionDataSource(lambda x: x**3 - 2 * x + 1, 0.0, 4.0, 7)

    rows = calculate_integration_comparison(source, 52.0, [2 / 3, 1 / 3])

    assert len(rows) == 6
    assert {row["method"] for row in rows} == {
        "CompositeTrapezoidal",
        "CompositeSimpsonOneThird",
        "CompositeSimpsonThreeEighths",
    }
    assert [row["n"] for row in rows] == [7, 7, 7, 13, 13, 13]
    assert [row["result"] for row in rows if row["method"] != "CompositeTrapezoidal"] == pytest.approx(
        [52, 52, 52, 52]
    )


def test_save_integration_result_writes_output_files(tmp_path, monkeypatch):
    source = FunctionDataSource(lambda x: x**2, 0.0, 2.0, 3)
    result = CompositeSimpsonOneThird(source).result_row(exact=8 / 3)
    monkeypatch.setattr("main.PROJECT_DIR", str(tmp_path))

    output_paths = save_integration_result(result, source)

    assert [path.rsplit("/", 1)[-1] for path in output_paths] == [
        "simpson_1_3_table.txt",
        "simpson_1_3_results.csv",
        "simpson_1_3_analysis.txt",
        "simpson_1_3_graph.png",
    ]
    assert "CompositeSimpsonOneThird" in (
        tmp_path / "output" / "simpson_1_3_table.txt"
    ).read_text(encoding="utf-8")
    assert (tmp_path / "output" / "simpson_1_3_graph.png").read_bytes().startswith(
        b"\x89PNG\r\n\x1a\n"
    )


@pytest.mark.parametrize(
    ("method_type", "points", "expected_prefix"),
    [
        (CompositeTrapezoidal, 5, "trapezoidal"),
        (CompositeSimpsonOneThird, 5, "simpson_1_3"),
        (CompositeSimpsonThreeEighths, 7, "simpson_3_8"),
    ],
)
def test_each_integration_rule_has_distinct_output_names(
    tmp_path, monkeypatch, method_type, points, expected_prefix
):
    source = FunctionDataSource(lambda x: x**3, 0.0, 3.0, points)
    result = method_type(source).result_row(exact=81 / 4)
    monkeypatch.setattr("main.PROJECT_DIR", str(tmp_path))

    output_paths = save_integration_result(result, source)

    assert all(expected_prefix in path for path in output_paths)
    assert all((tmp_path / "output" / path.rsplit("/", 1)[-1]).is_file() for path in output_paths)


def test_save_integration_comparison_uses_log_log_axes(tmp_path, monkeypatch):
    source = FunctionDataSource(lambda x: x**3, 0.0, 3.0, 7)
    rows = []
    for h in (0.5, 0.25):
        for method_type in (
            CompositeTrapezoidal,
            CompositeSimpsonOneThird,
            CompositeSimpsonThreeEighths,
        ):
            result = method_type(source).result_row(exact=81 / 4)
            rows.append(
                {
                    "h": h,
                    "method": result.method,
                    "n": result.n,
                    "result": result.result,
                    "exact": result.exact,
                    "absolute_error": result.absolute_error,
                }
            )
    monkeypatch.setattr(main, "PROJECT_DIR", str(tmp_path))

    original_subplots = plt.subplots
    captured = {}

    def capture_subplots(*args, **kwargs):
        figure, axes = original_subplots(*args, **kwargs)
        captured["axes"] = axes
        return figure, axes

    monkeypatch.setattr(plt, "subplots", capture_subplots)
    save_integration_comparison(rows, source)

    axes = captured["axes"].flatten()
    assert all(axis.get_xscale() == "log" for axis in axes)
    assert all(axis.get_yscale() == "log" for axis in axes)


def test_save_integration_comparison_writes_report_and_graph(tmp_path, monkeypatch):
    source = FunctionDataSource(lambda x: x**3, 0.0, 3.0, 7)
    rows = []
    for method_type in (
        CompositeTrapezoidal,
        CompositeSimpsonOneThird,
        CompositeSimpsonThreeEighths,
    ):
        result = method_type(source).result_row(exact=81 / 4)
        rows.append(
            {
                "h": 0.5,
                "method": result.method,
                "n": result.n,
                "result": result.result,
                "exact": result.exact,
                "absolute_error": result.absolute_error,
            }
        )
    monkeypatch.setattr(main, "PROJECT_DIR", str(tmp_path))

    output_paths = save_integration_comparison(rows, source)

    output_dir = tmp_path / "output"
    assert [path.rsplit("/", 1)[-1] for path in output_paths] == [
        "integration_comparison_table.txt",
        "integration_comparison_results.csv",
        "integration_comparison_analysis.txt",
        "integration_comparison_graph.png",
    ]
    assert "CompositeSimpsonOneThird" in (output_dir / "integration_comparison_table.txt").read_text(encoding="utf-8")
    assert (output_dir / "integration_comparison_graph.png").read_bytes().startswith(
        b"\x89PNG\r\n\x1a\n"
    )