"""Reconstruct Roodman's three introductory GWP charts through 2025.

The layout deliberately follows the published figures: a wide white canvas,
centered serif title, teal observations, and only sparse in-chart fit labels.
Red is reserved for the added 2020--2025 observations. This is a visual data
extension, not a rerun of the paper's stochastic estimation.
"""

from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np

from plot_gwp_extension import extension, original_gwp


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reference-output-sample" / "article-chart-reconstructions"

# Matches the published PNG aspect ratio (1323 × 789) and its sparse canvas.
WIDTH, HEIGHT = 1323, 789
LEFT, RIGHT, TOP, BOTTOM = 25, 8, 95, 20
PLOT_W, PLOT_H = WIDTH - LEFT - RIGHT, HEIGHT - TOP - BOTTOM
TEAL, RED, INK, PAPER = "#309ba6", "#cf3d32", "#353535", "#ffffff"


def line_path(xs: np.ndarray, ys: np.ndarray, sx, sy) -> str:
    return " ".join(
        ("M" if index == 0 else "L") + f" {sx(float(x)):.2f} {sy(float(y)):.2f}"
        for index, (x, y) in enumerate(zip(xs, ys))
    )


def chart(
    name: str,
    data_x: np.ndarray,
    data_y: np.ndarray,
    x_domain: tuple[float, float],
    y_domain: tuple[float, float],
    *,
    x_log: bool = False,
    y_log: bool = False,
    reverse_x: bool = False,
    fits: list[tuple[np.ndarray, np.ndarray, str, float, float]] | None = None,
    update_x: np.ndarray | None = None,
    update_y: np.ndarray | None = None,
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

    points = "".join(
        f'<circle cx="{sx(float(x)):.1f}" cy="{sy(float(y)):.1f}" r="4.0" class="history-point"/>'
        for x, y in zip(data_x, data_y)
    )
    series_path = f'<path d="{line_path(data_x, data_y, sx, sy)}" class="history-series"/>'
    fit_paths = "" if fits is None else "".join(
        f'<path d="{line_path(fit_x, fit_y, sx, sy)}" class="fit"/>'
        for fit_x, fit_y, _, _, _ in fits
    )
    fit_labels = "" if fits is None else "".join(
        f'<text x="{sx(label_x):.1f}" y="{sy(label_y):.1f}" class="fit-label">{escape(label)}</text>'
        for _, _, label, label_x, label_y in fits
    )
    update_path = "" if update_x is None or update_y is None else (
        f'<path d="{line_path(update_x, update_y, sx, sy)}" class="update-series"/>'
    )
    update_points = "" if update_x is None or update_y is None else "".join(
        f'<circle cx="{sx(float(x)):.1f}" cy="{sy(float(y)):.1f}" r="5.0" class="update-point"/>'
        for x, y in zip(update_x, update_y)
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">
<style>
  text {{ font-family: Georgia, "Times New Roman", serif; fill: {INK}; }}
  .title {{ font-size: 36px; }} .fit-label {{ font-size: 29px; }}
  .history-point {{ fill: {TEAL}; }}
  .history-series {{ fill: none; stroke: {TEAL}; stroke-width: 2.5; }}
  .update-series {{ fill: none; stroke: {RED}; stroke-width: 2.5; }}
  .update-point {{ fill: {RED}; stroke: {PAPER}; stroke-width: 1.25; }}
  .fit {{ fill: none; stroke: {TEAL}; stroke-width: 2.8; stroke-dasharray: 1 8; stroke-linecap: round; }}
</style>
<defs><clipPath id="plot-area"><rect x="{LEFT}" y="{TOP}" width="{PLOT_W}" height="{PLOT_H}"/></clipPath></defs>
<rect width="{WIDTH}" height="{HEIGHT}" fill="{PAPER}"/>
<text x="{WIDTH / 2:.1f}" y="52" text-anchor="middle" class="title">Gross world product, 10,000 BCE–2025 CE</text>
<g clip-path="url(#plot-area)">{fit_paths}{series_path}{points}{update_path}{update_points}</g>{fit_labels}
</svg>'''
    (OUT / name).write_text(svg)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    series = original_gwp()
    years = series["year"].to_numpy()
    gwp = series["gwp"].to_numpy()
    updated = extension(series)
    recent = updated.loc[updated["year"] >= 2020]
    recent_years = recent["year"].to_numpy()
    recent_gwp = recent["chain_linked_gwp"].to_numpy()

    exponential = np.polyfit(years, np.log(gwp), 1)
    fit_years = np.linspace(years.min(), years.max(), 600)
    exponential_fit = np.exp(np.polyval(exponential, fit_years))

    takeoff = 2047
    years_to_takeoff = takeoff - years
    power = np.polyfit(np.log(years_to_takeoff), np.log(gwp), 1)
    fit_distance = np.geomspace(years_to_takeoff.min(), years_to_takeoff.max(), 600)
    power_fit = np.exp(np.polyval(power, np.log(fit_distance)))

    # Publication-era reconstructions.
    chart("01-ordinary-axes.svg", years, gwp, (-10000, 2019), (0, 80000))
    chart(
        "02-log-y-exponential-fit.svg", years, gwp, (-10000, 2019), (1, 100000),
        y_log=True,
        fits=[(fit_years, exponential_fit, "Exponential fit", -2000, 13000)],
    )
    chart(
        "03-transformed-time-power-fit.svg", years_to_takeoff, gwp, (20, 15000), (1, 100000),
        x_log=True, y_log=True, reverse_x=True,
        fits=[
            (2047 - fit_years, exponential_fit, "Exponential fit", 3500, 210),
            (fit_distance, power_fit, "Power law fit", 120, 14000),
        ],
    )

    # The corresponding visual extensions: no new framing, just red observations.
    chart(
        "01-ordinary-axes-through-2025.svg", years, gwp, (-10000, 2025), (0, 100000),
        update_x=recent_years, update_y=recent_gwp,
    )
    chart(
        "02-log-y-exponential-fit-through-2025.svg", years, gwp, (-10000, 2025), (1, 100000),
        y_log=True,
        fits=[(fit_years, exponential_fit, "Exponential fit", -2000, 13000)],
        update_x=recent_years, update_y=recent_gwp,
    )
    chart(
        "03-transformed-time-power-fit-through-2025.svg", years_to_takeoff, gwp, (20, 15000), (1, 100000),
        x_log=True, y_log=True, reverse_x=True,
        fits=[
            (2047 - fit_years, exponential_fit, "Exponential fit", 3500, 210),
            (fit_distance, power_fit, "Power law fit", 120, 14000),
        ],
        update_x=takeoff - recent_years, update_y=recent_gwp,
    )


if __name__ == "__main__":
    main()
