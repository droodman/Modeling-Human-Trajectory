"""Reconstruct Roodman's three introductory GWP charts through 2025.

The layout deliberately follows the published figures: a wide white canvas,
centered serif title, teal observations, and only sparse in-chart fit labels.
Red is reserved for the added 2020--2025 observations. This is a visual data
extension, not a rerun of the paper's stochastic estimation.
"""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "GWP.xlsx"
WORLD_BANK = ROOT / "data-update" / "world-bank-gwp-2019-2025.json"
OUT = ROOT / "reference-output-sample"

# Matches the published PNG aspect ratio (1323 × 789) and its sparse canvas.
WIDTH, HEIGHT = 1323, 789
LEFT, RIGHT, TOP, BOTTOM = 25, 8, 95, 20
PLOT_W, PLOT_H = WIDTH - LEFT - RIGHT, HEIGHT - TOP - BOTTOM
TEAL, RED, INK, PAPER = "#309ba6", "#cf3d32", "#353535", "#ffffff"


def original_gwp() -> pd.DataFrame:
    """Recreate the GWP series used by the first PrepData pass in Model GWP.do."""
    data = pd.read_excel(WORKBOOK, sheet_name="Data", header=1)
    frame = pd.DataFrame(
        {
            "year": data["Year"].astype(float),
            "maddison_population": data["Pop, Maddison 2010"],
            "mcevedy_population": data["Pop, McEvedy & Jones"],
            "un_population": data["Pop, UN"],
            "gwp_per_capita": data["GWP/cap, Maddison 2010"],
            "weo_growth": data["GWP growth, WEO"],
        }
    )
    frame.loc[frame["year"] == 1, "year"] = 0
    frame = frame.loc[frame["year"] >= -10000].copy().reset_index(drop=True)
    frame["population"] = np.where(
        frame["year"] < 1500, frame["mcevedy_population"], frame["maddison_population"]
    )
    last_population = frame["population"].last_valid_index()
    base_population = frame.loc[last_population, "population"]
    base_un_population = frame.loc[last_population, "un_population"]
    after_population = frame.index > last_population
    frame.loc[after_population, "population"] = (
        base_population * frame.loc[after_population, "un_population"] / base_un_population
    )
    for start, end in ((1000, 1500), (0, 1000), (-10000, 0)):
        start_row = frame.index[frame["year"] == start][0]
        end_row = frame.index[frame["year"] == end][0]
        start_value = 400.0 if start == -10000 else frame.loc[start_row, "gwp_per_capita"]
        end_value = frame.loc[end_row, "gwp_per_capita"]
        mask = frame["year"].between(start, end)
        share = (
            np.log(frame.loc[mask, "population"]) - np.log(frame.loc[start_row, "population"])
        ) / (np.log(frame.loc[end_row, "population"]) - np.log(frame.loc[start_row, "population"]))
        frame.loc[mask, "gwp_per_capita"] = start_value * (end_value / start_value) ** share
    frame["gwp"] = frame["gwp_per_capita"] * frame["population"] / 1000
    for index in frame.index:
        if frame.loc[index, "year"] >= 2000 and pd.isna(frame.loc[index, "gwp"]):
            frame.loc[index, "gwp"] = frame.loc[index - 1, "gwp"] * (1 + frame.loc[index, "weo_growth"])
    return frame.loc[frame["gwp"].notna(), ["year", "gwp"]]


def extension(history: pd.DataFrame) -> pd.DataFrame:
    """Chain-link World Bank PPP-GDP growth to the workbook's 2019 GWP level."""
    payload = json.loads(WORLD_BANK.read_text())
    records = payload["response"][1]
    world_bank = pd.DataFrame(
        {
            "year": [int(row["date"]) for row in records],
            "world_bank_ppp_gdp": [row["value"] for row in records],
        }
    ).sort_values("year")
    wb_2019 = world_bank.loc[world_bank["year"] == 2019, "world_bank_ppp_gdp"].iloc[0]
    original_2019 = history.loc[history["year"] == 2019, "gwp"].iloc[0]
    world_bank["chain_linked_gwp"] = original_2019 * world_bank["world_bank_ppp_gdp"] / wb_2019
    return world_bank


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
    # Keep the coefficients fitted to the original through-2019 observations,
    # but project the unchanged fit through the added 2025 endpoint.
    fit_years = np.linspace(years.min(), recent_years.max(), 600)
    exponential_fit = np.exp(np.polyval(exponential, fit_years))

    takeoff = 2047
    years_to_takeoff = takeoff - years
    power = np.polyfit(np.log(years_to_takeoff), np.log(gwp), 1)
    fit_distance = np.geomspace(
        takeoff - recent_years.max(), years_to_takeoff.max(), 600
    )
    power_fit = np.exp(np.polyval(power, np.log(fit_distance)))

    # Three visual extensions in the published layout: no new framing, just red observations.
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
