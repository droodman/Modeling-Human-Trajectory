"""Reconstruct the three introductory GWP charts from Roodman's article.

The series comes from the data-preparation logic copied from `Model GWP.do` in
`plot_gwp_extension.py`; this script changes axes and fits reference lines only.
It does not rerun the paper's stochastic estimation.
"""

from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np

from plot_gwp_extension import original_gwp


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reference-output-sample" / "article-chart-reconstructions"
WIDTH, HEIGHT = 1200, 1200
LEFT, RIGHT, TOP, BOTTOM = 140, 65, 115, 200
PLOT_W, PLOT_H = WIDTH - LEFT - RIGHT, HEIGHT - TOP - BOTTOM
RED, INK, GRID, PAPER = "#e53b2c", "#222222", "#dedbd5", "#fbfaf7"


def line_path(xs: np.ndarray, ys: np.ndarray, sx, sy) -> str:
    return " ".join(
        ("M" if index == 0 else "L") + f" {sx(float(x)):.2f} {sy(float(y)):.2f}"
        for index, (x, y) in enumerate(zip(xs, ys))
    )


def chart(
    name: str,
    title: str,
    subtitle: str,
    data_x: np.ndarray,
    data_y: np.ndarray,
    x_domain: tuple[float, float],
    y_domain: tuple[float, float],
    x_ticks: list[tuple[float, str]],
    y_ticks: list[tuple[float, str]],
    x_label: str,
    y_label: str,
    *,
    x_log: bool = False,
    y_log: bool = False,
    reverse_x: bool = False,
    fit: tuple[np.ndarray, np.ndarray] | None = None,
) -> None:
    def scale(value: float, domain: tuple[float, float], span: float, log: bool) -> float:
        lo, hi = domain
        if log:
            value, lo, hi = np.log(value), np.log(lo), np.log(hi)
        return span * (value - lo) / (hi - lo)

    def sx(value: float) -> float:
        fraction = scale(value, x_domain, PLOT_W, x_log)
        if reverse_x:
            fraction = PLOT_W - fraction
        return LEFT + fraction

    def sy(value: float) -> float:
        return TOP + PLOT_H - scale(value, y_domain, PLOT_H, y_log)

    grid_lines, labels = [], []
    for value, label in x_ticks:
        x = sx(value)
        grid_lines.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{TOP + PLOT_H}" class="grid"/>')
        labels.append(f'<text x="{x:.1f}" y="{TOP + PLOT_H + 31}" text-anchor="middle" class="tick">{escape(label)}</text>')
    for value, label in y_ticks:
        y = sy(value)
        grid_lines.append(f'<line x1="{LEFT}" y1="{y:.1f}" x2="{LEFT + PLOT_W}" y2="{y:.1f}" class="grid"/>')
        labels.append(f'<text x="{LEFT - 15}" y="{y + 5:.1f}" text-anchor="end" class="tick">{escape(label)}</text>')

    points = "".join(
        f'<circle cx="{sx(float(x)):.1f}" cy="{sy(float(y)):.1f}" r="3.0" class="point"/>'
        for x, y in zip(data_x, data_y)
    )
    series_path = f'<path d="{line_path(data_x, data_y, sx, sy)}" class="series"/>'
    fit_path = "" if fit is None else f'<path d="{line_path(fit[0], fit[1], sx, sy)}" class="fit"/>'
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">
<style>
  text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill: {INK}; }}
  .title {{ font-size: 23px; font-weight: 700; }} .subtitle {{ font-size: 13px; fill: #5b5b5b; }}
  .tick {{ font-size: 14px; fill: #555; }} .grid {{ stroke: {GRID}; stroke-width: 1; }}
  .axis {{ stroke: {INK}; stroke-width: 1.5; }} .point {{ fill: {RED}; }}
  .series {{ fill: none; stroke: {RED}; stroke-width: 1.4; stroke-opacity: 0.55; }}
  .fit {{ fill: none; stroke: {RED}; stroke-width: 2.2; stroke-dasharray: 3 7; }} .note {{ font-size: 11px; fill: #5b5b5b; }}
</style>
<defs><clipPath id="plot-area"><rect x="{LEFT}" y="{TOP}" width="{PLOT_W}" height="{PLOT_H}"/></clipPath></defs>
<rect width="100%" height="100%" fill="{PAPER}"/>
<text x="{LEFT}" y="45" class="title">{escape(title)}</text>
<text x="{LEFT}" y="73" class="subtitle">{escape(subtitle)}</text>
{''.join(grid_lines)}
<line x1="{LEFT}" y1="{TOP}" x2="{LEFT}" y2="{TOP + PLOT_H}" class="axis"/>
<line x1="{LEFT}" y1="{TOP + PLOT_H}" x2="{LEFT + PLOT_W}" y2="{TOP + PLOT_H}" class="axis"/>
<g clip-path="url(#plot-area)">{series_path}{fit_path}{points}</g>{''.join(labels)}
<text x="{LEFT + PLOT_W / 2:.1f}" y="{TOP + PLOT_H + 77}" text-anchor="middle" class="tick">{escape(x_label)}</text>
<text x="33" y="{TOP + PLOT_H / 2:.1f}" transform="rotate(-90 33 {TOP + PLOT_H / 2:.1f})" text-anchor="middle" class="tick">{escape(y_label)}</text>
<text x="{LEFT}" y="{HEIGHT - 48}" class="note">Reconstructed from Roodman GWP.xlsx and the data-preparation logic in Model GWP.do.</text>
<text x="{LEFT}" y="{HEIGHT - 27}" class="note">The dotted reference line is a simple least-squares fit in the chart’s displayed coordinates.</text>
</svg>'''
    (OUT / name).write_text(svg)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    series = original_gwp()
    years = series["year"].to_numpy()
    gwp = series["gwp"].to_numpy()

    chart(
        "01-ordinary-axes.svg",
        "Global world product, 10,000 BCE–2019",
        "The same reconstructed series on ordinary linear axes: the familiar hockey stick.",
        years,
        gwp,
        (-10500, 2500),
        (0, 80000),
        [(-10000, "10,000 BCE"), (-5000, "5,000 BCE"), (0, "1 CE"), (1000, "1000"), (2000, "2000")],
        [(0, "0"), (20000, "20,000"), (40000, "40,000"), (60000, "60,000"), (80000, "80,000")],
        "Year",
        "GWP (1990 international-$ billions)",
    )

    exponential = np.polyfit(years, np.log(gwp), 1)
    fit_years = np.linspace(years.min(), years.max(), 600)
    chart(
        "02-log-y-exponential-fit.svg",
        "The same series on a logarithmic vertical axis",
        "Equal vertical spacing means a tenfold increase; the dotted line is the best exponential fit.",
        years,
        gwp,
        (-10500, 2500),
        (1, 100000),
        [(-10000, "10,000 BCE"), (-5000, "5,000 BCE"), (0, "1 CE"), (1000, "1000"), (2000, "2000")],
        [(1, "$1b"), (10, "$10b"), (100, "$100b"), (1000, "$1t"), (10000, "$10t"), (100000, "$100t")],
        "Year",
        "GWP (1990 international dollars, log scale)",
        y_log=True,
        fit=(fit_years, np.exp(np.polyval(exponential, fit_years))),
    )

    takeoff = 2047
    years_to_takeoff = takeoff - years
    power = np.polyfit(np.log(years_to_takeoff), np.log(gwp), 1)
    fit_distance = np.geomspace(years_to_takeoff.min(), years_to_takeoff.max(), 600)
    chart(
        "03-transformed-time-power-fit.svg",
        "The series on Roodman’s transformed time axis",
        "Equal horizontal spacing means a tenfold reduction in years until 2047; the dotted line is the best power-law fit.",
        years_to_takeoff,
        gwp,
        (20, 15000),
        (1, 100000),
        [(10000, "10,000"), (1000, "1,000"), (100, "100"), (10, "10")],
        [(1, "$1b"), (10, "$10b"), (100, "$100b"), (1000, "$1t"), (10000, "$10t"), (100000, "$100t")],
        "Years until 2047 (log scale, reversed)",
        "GWP (1990 international dollars, log scale)",
        x_log=True,
        y_log=True,
        reverse_x=True,
        fit=(fit_distance, np.exp(np.polyval(power, np.log(fit_distance)))),
    )


if __name__ == "__main__":
    main()
