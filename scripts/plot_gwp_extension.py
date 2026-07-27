"""Extend Roodman's transformed GWP chart beyond the workbook's 2019 endpoint.

The historical side mirrors the core data-preparation logic in `PrepData` from
`Model GWP.do` for the GWP chart. The update is *not* a fresh re-estimation:
World Bank global PPP GDP levels are chain-linked to the model's 2019 endpoint,
so its post-2019 growth is visible in the same visual scale without claiming
that its price base or source construction is identical to the original model.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "GWP.xlsx"
WORLD_BANK = ROOT / "data-update" / "world-bank-world-gdp-ppp-constant-2021-2019-2025.json"
CSV_OUT = ROOT / "data-update" / "gwp-extension-2019-2025.csv"
PLOT_OUT = ROOT / "reference-output-sample" / "gwp-through-2025-chain-linked.svg"


def original_gwp() -> pd.DataFrame:
    """Recreate the GWP series used by the first `PrepData` call in the do-file."""
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
    # Stata changes 1 CE to year 0 before calculating the transformed axis.
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

    # The three pre-1500 segments interpolate GWP per person in log population.
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
    # This is Stata's forward-fill of GWP using the then-current WEO growth rate.
    for idx in frame.index:
        if frame.loc[idx, "year"] >= 2000 and pd.isna(frame.loc[idx, "gwp"]):
            frame.loc[idx, "gwp"] = frame.loc[idx - 1, "gwp"] * (1 + frame.loc[idx, "weo_growth"])
    return frame.loc[frame["gwp"].notna(), ["year", "gwp"]]


def extension(history: pd.DataFrame) -> pd.DataFrame:
    """Chain-link World Bank PPP GDP levels to the original 2019 GWP endpoint."""
    payload = json.loads(WORLD_BANK.read_text())
    records = payload[1]
    wb = pd.DataFrame(
        {"year": [int(row["date"]) for row in records], "world_bank_ppp_gdp": [row["value"] for row in records]}
    ).sort_values("year")
    wb = wb.loc[wb["world_bank_ppp_gdp"].notna()].copy()
    wb_2019 = wb.loc[wb["year"] == 2019, "world_bank_ppp_gdp"].iloc[0]
    original_2019 = history.loc[history["year"] == 2019, "gwp"].iloc[0]
    wb["chain_linked_gwp"] = original_2019 * wb["world_bank_ppp_gdp"] / wb_2019
    wb["series"] = np.where(wb["year"] == 2019, "join point", "World Bank update")
    return wb


def draw(history: pd.DataFrame, update: pd.DataFrame) -> None:
    """Write a self-contained SVG to avoid a chart-library dependency."""
    width, height = 1440, 860
    left, right, top, bottom = 145, 80, 105, 155
    plot_w, plot_h = width - left - right, height - top - bottom
    x_lo, x_hi = 20, 15000
    y_lo, y_hi = 1, 100000

    def sx(year: float) -> float:
        distance = 2050 - year
        return left + plot_w * (np.log(x_hi) - np.log(distance)) / (np.log(x_hi) - np.log(x_lo))

    def sy(value: float) -> float:
        return top + plot_h * (np.log(y_hi) - np.log(value)) / (np.log(y_hi) - np.log(y_lo))

    def path(rows: pd.DataFrame, value_col: str) -> str:
        return " ".join(
            ("M" if position == 0 else "L") + f" {sx(row.year):.1f} {sy(getattr(row, value_col)):.1f}"
            for position, row in enumerate(rows.itertuples(index=False))
        )

    historical_path = path(history, "gwp")
    new = update.loc[update["year"] >= 2020]
    update_path = path(new, "chain_linked_gwp")
    major_x = (10000, 1000, 100, 20)
    major_y = (1, 10, 100, 1000, 10000, 100000)
    grid = []
    labels = []
    for value in major_x:
        x = sx(2050 - value)
        grid.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" class="grid"/>')
        labels.append(f'<text x="{x:.1f}" y="{top + plot_h + 28}" text-anchor="middle" class="tick">{value:,}</text>')
    for value in major_y:
        y = sy(value)
        grid.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>')
        labels.append(f'<text x="{left - 14}" y="{y + 5:.1f}" text-anchor="end" class="tick">{value:,}</text>')
    historical_dots = "".join(
        f'<circle cx="{sx(row.year):.1f}" cy="{sy(row.gwp):.1f}" r="2.2" class="history-dot"/>'
        for row in history.itertuples(index=False)
    )
    update_dots = "".join(
        f'<circle cx="{sx(row.year):.1f}" cy="{sy(row.chain_linked_gwp):.1f}" r="5" class="update-dot"/>'
        for row in new.itertuples(index=False)
    )
    join_value = update.loc[update["year"] == 2019, "chain_linked_gwp"].iloc[0]
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
  text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill: #252525; }}
  .title {{ font-size: 25px; font-weight: 700; }} .subtitle {{ font-size: 14px; fill: #5b5b5b; }}
  .tick {{ font-size: 14px; fill: #555; }} .grid {{ stroke: #dedbd5; stroke-width: 1; }}
  .axis {{ stroke: #282828; stroke-width: 1.4; }} .history {{ fill: none; stroke: #1f77b4; stroke-width: 1.6; }}
  .history-dot {{ fill: #1f77b4; }} .update {{ fill: none; stroke: #d63d31; stroke-width: 2.8; }}
  .update-dot {{ fill: #d63d31; stroke: #fbfaf7; stroke-width: 1.4; }} .legend {{ font-size: 13px; }} .note {{ font-size: 12px; fill: #5b5b5b; }}
</style>
<rect width="100%" height="100%" fill="#fbfaf7"/>
<text x="{left}" y="43" class="title">Roodman GWP series extended with post-2019 global GDP data</text>
<text x="{left}" y="68" class="subtitle">Blue: original workbook reconstruction through 2019. Red: World Bank global PPP GDP, rescaled to meet it at 2019.</text>
{''.join(grid)}
<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>
<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" class="axis"/>
<path d="{historical_path}" class="history"/>{historical_dots}
<path d="{update_path}" class="update"/>{update_dots}
<circle cx="{sx(2019):.1f}" cy="{sy(join_value):.1f}" r="5" fill="#fbfaf7" stroke="#d63d31" stroke-width="2"/>
{''.join(labels)}
<text x="{left + plot_w / 2:.1f}" y="{top + plot_h + 72}" text-anchor="middle" class="tick">Years until 2050 (log scale, reversed)</text>
<text x="31" y="{top + plot_h / 2:.1f}" transform="rotate(-90 31 {top + plot_h / 2:.1f})" text-anchor="middle" class="tick">GWP (original model scale; $ billions, log scale)</text>
<line x1="{left + 15}" y1="{top + 20}" x2="{left + 49}" y2="{top + 20}" class="history"/><circle cx="{left + 32}" cy="{top + 20}" r="2.4" class="history-dot"/>
<text x="{left + 58}" y="{top + 25}" class="legend">Original reconstructed series through 2019</text>
<line x1="{left + 15}" y1="{top + 44}" x2="{left + 49}" y2="{top + 44}" class="update"/><circle cx="{left + 32}" cy="{top + 44}" r="4" class="update-dot"/>
<text x="{left + 58}" y="{top + 49}" class="legend">2020–2025 World Bank update, chain-linked at 2019</text>
<text x="{left}" y="{height - 45}" class="note">Sources: Roodman GWP.xlsx (May 2020 archive); World Bank WDI NY.GDP.MKTP.PP.KD, retrieved 2026-07-27.</text>
<text x="{left}" y="{height - 25}" class="note">The red extension is a visual update, not a rerun or re-estimation of Roodman’s stochastic model.</text>
</svg>'''
    PLOT_OUT.write_text(svg)


if __name__ == "__main__":
    historical = original_gwp()
    updated = extension(historical)
    CSV_OUT.parent.mkdir(exist_ok=True)
    updated.to_csv(CSV_OUT, index=False)
    draw(historical, updated)
    print(updated[["year", "chain_linked_gwp", "series"]].to_string(index=False))
