# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy>=2.0",
#   "openpyxl>=3.1",
#   "pandas>=2.2",
#   "scipy>=1.14",
# ]
# ///
"""Extend Roodman's GWP prediction-percentile analysis through 2025.

The script ports the univariate Feller-diffusion likelihood and simulation
used by Model GWP.do/asdf. It validates the port against archived full-sample
estimates and the published 2019 rolling forecast before calculating rolling
ten-year forecasts and two alternative forecast constructions.
"""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special


ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "GWP.xlsx"
WORLD_BANK = ROOT / "data-update" / "world-bank-gwp-2019-2025.json"
MAIN_OUTPUT = (
    ROOT
    / "reference-output-sample"
    / "gwp-prediction-percentiles-through-2025.svg"
)
ANNUAL_ALTERNATIVE_OUTPUT = (
    ROOT
    / "reference-output-sample"
    / "gwp-prediction-percentiles-alternative-annual.svg"
)
FIXED_ALTERNATIVE_OUTPUT = (
    ROOT
    / "reference-output-sample"
    / "gwp-prediction-percentiles-alternative-fixed-2010.svg"
)

# Digitized from the author's published 1440 x 1047 PNG. Keeping these values
# fixed prevents a fresh Monte Carlo draw from visually moving the old dots.
PUBLISHED_PERCENTILES = {
    1600: 0.8570,
    1700: 0.6576,
    1820: 0.9404,
    1870: 0.9361,
    1913: 0.9816,
    1940: 0.7822,
    1950: 0.2524,
    1960: 0.8938,
    1970: 0.7671,
    1980: 0.4399,
    1990: 0.2589,
    2000: 0.2243,
    2010: 0.2568,
    2019: 0.2080,
}

EXPECTED_THETA = np.array([-12.661317, 1.857147e-5, -23.778794, -1.812893])
EXPECTED_STANDARD_ERRORS = np.array([0.280591, 6.87266e-5, 7.43875, 0.161617])


def original_gwp() -> pd.DataFrame:
    """Recreate the first PrepData pass in Model GWP.do."""
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
        frame["year"] < 1500,
        frame["mcevedy_population"],
        frame["maddison_population"],
    )
    last_population = frame["population"].last_valid_index()
    base_population = frame.loc[last_population, "population"]
    base_un_population = frame.loc[last_population, "un_population"]
    after_population = frame.index > last_population
    frame.loc[after_population, "population"] = (
        base_population
        * frame.loc[after_population, "un_population"]
        / base_un_population
    )

    for start, end in ((1000, 1500), (0, 1000), (-10000, 0)):
        start_row = frame.index[frame["year"] == start][0]
        end_row = frame.index[frame["year"] == end][0]
        start_value = (
            400.0 if start == -10000 else frame.loc[start_row, "gwp_per_capita"]
        )
        end_value = frame.loc[end_row, "gwp_per_capita"]
        mask = frame["year"].between(start, end)
        share = (
            np.log(frame.loc[mask, "population"])
            - np.log(frame.loc[start_row, "population"])
        ) / (
            np.log(frame.loc[end_row, "population"])
            - np.log(frame.loc[start_row, "population"])
        )
        frame.loc[mask, "gwp_per_capita"] = start_value * (
            end_value / start_value
        ) ** share

    frame["gwp"] = frame["gwp_per_capita"] * frame["population"] / 1000
    for index in frame.index:
        if frame.loc[index, "year"] >= 2000 and pd.isna(frame.loc[index, "gwp"]):
            frame.loc[index, "gwp"] = frame.loc[index - 1, "gwp"] * (
                1 + frame.loc[index, "weo_growth"]
            )
    return frame.loc[frame["gwp"].notna(), ["year", "gwp"]]


def updated_gwp(history: pd.DataFrame) -> pd.DataFrame:
    """Chain-link World Bank PPP-GDP growth to the reconstructed 2019 level."""
    payload = json.loads(WORLD_BANK.read_text())
    records = payload["response"][1]
    world_bank = pd.DataFrame(
        {
            "year": [int(row["date"]) for row in records],
            "world_bank_ppp_gdp": [row["value"] for row in records],
        }
    ).sort_values("year")
    wb_2019 = world_bank.loc[
        world_bank["year"] == 2019, "world_bank_ppp_gdp"
    ].iloc[0]
    original_2019 = history.loc[history["year"] == 2019, "gwp"].iloc[0]
    world_bank["gwp"] = (
        original_2019 * world_bank["world_bank_ppp_gdp"] / wb_2019
    )
    return world_bank[["year", "gwp"]]


def model_frame(include_update: bool = False) -> pd.DataFrame:
    history = original_gwp().rename(columns={"gwp": "GWP"})
    history = history.loc[
        (history["year"] >= -10000)
        & (
            (history["year"] <= 1950)
            | history["year"].isin(
                [1960, 1970, 1980, 1990, 2000, 2010, 2019]
            )
        )
    ].copy()
    if include_update:
        recent = updated_gwp(original_gwp())
        recent = recent.loc[recent["year"].between(2020, 2025)].rename(
            columns={"gwp": "GWP"}
        )
        history = pd.concat([history, recent], ignore_index=True)

    # HYDE uncertainty knots from PrepData in Model GWP.do.
    knots_year = np.array([-10000.0, 0.0, 1700.0, 1900.0, 2000.0, 2025.0])
    knots_sd = np.array([1.0, 0.75, 0.25, 0.05, 0.01, 0.01])
    history["HYDEsd"] = np.interp(history["year"], knots_year, knots_sd)
    history["weight"] = 1.0 / (1.0 + 2.0 * history["HYDEsd"] ** 2)
    return history.sort_values("year").reset_index(drop=True)


def unpack(q: np.ndarray) -> np.ndarray:
    """Undo optimizer scaling for [ln(a), b, nu, gamma]."""
    return np.array([q[0], q[1] * 1e-4, q[2], q[3]], dtype=float)


def transition_logpdf(
    y0: np.ndarray,
    y1: np.ndarray,
    dt: np.ndarray,
    theta: np.ndarray,
) -> np.ndarray:
    lna, b, nu, gamma = theta
    if not (
        -100.0 < lna < 20.0
        and -1000.0 < nu < 0.0
        and -50.0 < gamma < -1e-4
    ):
        return np.full_like(y0, -np.inf, dtype=float)

    a = np.exp(lna)
    bt = b * dt
    c_scale = (
        1.0 / (a * dt)
        if abs(b) < 1e-11
        else (b / np.expm1(bt)) / a
    )
    if np.any(~np.isfinite(c_scale)) or np.any(c_scale <= 0):
        return np.full_like(y0, -np.inf, dtype=float)

    log_x0 = np.log(y0) / gamma
    log_x = np.log(y1) / gamma
    log_c = np.log(c_scale)
    log_lam = log_x0 + log_c + bt
    log_z = log_x + log_c
    lam = np.exp(log_lam)
    z = np.exp(log_z)
    argument = 2.0 * np.sqrt(lam * z)
    scaled_bessel = special.ive(-nu, argument)
    log_bessel = np.log(scaled_bessel) + np.abs(argument)
    log_feller = (
        -lam
        - z
        + 0.5 * nu * (log_z - log_lam)
        + log_bessel
    )
    log_jacobian = (
        log_c
        - np.log(abs(gamma))
        + (1.0 / gamma - 1.0) * np.log(y1)
    )
    return log_feller + log_jacobian


def objective(q: np.ndarray, frame: pd.DataFrame) -> float:
    y = frame["GWP"].to_numpy(float)
    years = frame["year"].to_numpy(float)
    weights = frame["weight"].to_numpy(float)[1:]
    log_likelihood = transition_logpdf(
        y[:-1], y[1:], np.diff(years), unpack(q)
    )
    if np.any(~np.isfinite(log_likelihood)):
        return 1e100
    return -float(np.dot(weights, log_likelihood))


def fit(frame: pd.DataFrame, q0: np.ndarray | None = None):
    if q0 is None:
        q0 = np.array([-12.66, 0.186, -23.78, -1.813])
    starts = [
        q0,
        np.array([-13.0, 0.2, -20.0, -1.8]),
        np.array([-10.0, 0.0, -5.0, -1.0]),
    ]
    results = [
        optimize.minimize(
            objective,
            start,
            args=(frame,),
            method="Nelder-Mead",
            options={"maxiter": 10000, "xatol": 1e-9, "fatol": 1e-7},
        )
        for start in starts
    ]
    best = min(results, key=lambda result: result.fun)
    polished = optimize.minimize(
        objective,
        best.x,
        args=(frame,),
        method="L-BFGS-B",
        bounds=[(-30, 2), (-100, 100), (-300, -1e-5), (-20, -0.05)],
        options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-8, "maxls": 100},
    )
    return polished if polished.fun <= best.fun else best


def finite_hessian(function, point: np.ndarray, steps: np.ndarray) -> np.ndarray:
    size = len(point)
    hessian = np.empty((size, size), dtype=float)
    at_point = function(point)
    for i in range(size):
        ei = np.zeros(size)
        ei[i] = steps[i]
        hessian[i, i] = (
            function(point + ei) - 2.0 * at_point + function(point - ei)
        ) / steps[i] ** 2
        for j in range(i):
            ej = np.zeros(size)
            ej[j] = steps[j]
            value = (
                function(point + ei + ej)
                - function(point + ei - ej)
                - function(point - ei + ej)
                + function(point - ei - ej)
            ) / (4.0 * steps[i] * steps[j])
            hessian[i, j] = hessian[j, i] = value
    return hessian


def covariance_q(result, frame: pd.DataFrame) -> np.ndarray:
    mean_weight = frame["weight"].to_numpy(float)[1:].mean()
    raw = finite_hessian(
        lambda q: objective(q, frame),
        result.x,
        np.array([2e-4, 2e-4, 2e-3, 2e-4]),
    )
    # Stata normalizes analytic weights to have mean one.
    hessian = raw / mean_weight
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    clipped = np.maximum(eigenvalues, np.max(eigenvalues) * 1e-12)
    return (eigenvectors * (1.0 / clipped)) @ eigenvectors.T


def simulate_percentile(
    fit_result,
    frame: pd.DataFrame,
    start_value: float,
    observed_value: float,
    horizon: float,
    *,
    paths: int = 10000,
    steps: int = 10000,
    seed: int,
) -> float:
    covariance = covariance_q(fit_result, frame)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    square_root = eigenvectors @ np.diag(
        np.sqrt(np.maximum(eigenvalues, 0.0))
    )
    rng = np.random.default_rng(seed)
    draws_q = (
        fit_result.x[:, None]
        + square_root @ rng.normal(size=(4, paths))
    )
    draws = np.vstack(
        [draws_q[0], draws_q[1] * 1e-4, draws_q[2], draws_q[3]]
    )
    lna, b, nu, gamma = draws
    dt = horizon / steps
    a = np.exp(lna)
    drift_constant = dt * (1.0 + nu) * a
    noise_scale = np.sqrt(2.0 * dt * a)
    drift_factor = 1.0 + dt * b
    state = np.power(start_value, 1.0 / gamma)

    for _ in range(steps):
        active = state != 0.0
        state = (
            drift_factor * state
            + drift_constant * active
            + np.sqrt(state) * rng.normal(size=paths) * noise_scale
        )
        state[state < 0.0] = 0.0

    finite = state > 0.0
    final = np.full(paths, np.inf)
    final[finite] = np.power(state[finite], gamma[finite])
    return float(np.mean(final < observed_value))


def validate_port(frame: pd.DataFrame):
    full_fit = fit(frame)
    covariance = covariance_q(full_fit, frame)
    theta_scale = np.diag([1.0, 1e-4, 1.0, 1.0])
    standard_errors = np.sqrt(
        np.diag(theta_scale @ covariance @ theta_scale)
    )
    np.testing.assert_allclose(
        unpack(full_fit.x), EXPECTED_THETA, rtol=7e-4, atol=1e-7
    )
    np.testing.assert_allclose(
        standard_errors, EXPECTED_STANDARD_ERRORS, rtol=2e-3, atol=1e-8
    )

    prior_2019 = frame.loc[frame["year"] <= 2010].copy()
    fit_2019 = fit(prior_2019, full_fit.x)
    percentile_2019 = simulate_percentile(
        fit_2019,
        prior_2019,
        start_value=float(prior_2019.iloc[-1]["GWP"]),
        observed_value=float(
            frame.loc[frame["year"] == 2019, "GWP"].iloc[0]
        ),
        horizon=9.0,
        seed=902381476,
    )
    if abs(percentile_2019 - 0.2071) > 0.015:
        raise AssertionError(
            f"2019 holdout was {percentile_2019:.4f}; expected about 0.2071"
        )
    return full_fit, standard_errors, percentile_2019


def with_weights(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach the observation weights reconstructed in model_frame."""
    knots_year = np.array([-10000.0, 0.0, 1700.0, 1900.0, 2000.0, 2025.0])
    knots_sd = np.array([1.0, 0.75, 0.25, 0.05, 0.01, 0.01])
    weighted = frame.copy()
    weighted["HYDEsd"] = np.interp(weighted["year"], knots_year, knots_sd)
    weighted["weight"] = 1.0 / (1.0 + 2.0 * weighted["HYDEsd"] ** 2)
    return weighted.sort_values("year").reset_index(drop=True)


def annual_gwp() -> pd.DataFrame:
    """Combine workbook and chain-linked World Bank annual GWP levels."""
    workbook = original_gwp().rename(columns={"gwp": "GWP"})
    workbook = workbook.loc[
        workbook["year"].between(2011, 2019), ["year", "GWP"]
    ]
    world_bank = updated_gwp(original_gwp()).rename(columns={"gwp": "GWP"})
    world_bank = world_bank.loc[
        world_bank["year"].between(2020, 2025), ["year", "GWP"]
    ]
    return (
        pd.concat([workbook, world_bank], ignore_index=True)
        .drop_duplicates("year", keep="last")
        .sort_values("year")
        .reset_index(drop=True)
    )


def calculate_rolling_decade(q_start: np.ndarray) -> pd.DataFrame:
    """Calculate overlapping ten-year forecasts ending in 2020-2025."""
    base = model_frame(include_update=False)
    base = base.loc[base["year"] <= 2010, ["year", "GWP"]].copy()
    annual = annual_gwp()
    levels = pd.concat(
        [base.loc[base["year"] == 2010, ["year", "GWP"]], annual],
        ignore_index=True,
    ).set_index("year")["GWP"]

    missing = set(range(2010, 2026)).difference(int(year) for year in levels.index)
    if missing:
        raise RuntimeError(f"Missing annual GWP observations: {sorted(missing)}")

    rows = []
    for target_year in range(2020, 2026):
        origin_year = target_year - 10
        added = annual.loc[annual["year"].between(2011, origin_year)]
        prior = with_weights(
            pd.concat([base, added], ignore_index=True).drop_duplicates(
                "year", keep="last"
            )
        )
        fitted = fit(prior, q_start)
        percentile = simulate_percentile(
            fitted,
            prior,
            start_value=float(levels.loc[origin_year]),
            observed_value=float(levels.loc[target_year]),
            horizon=10.0,
            seed=902381476 + target_year,
        )
        rows.append(
            {
                "year": target_year,
                "fit_through": origin_year,
                "horizon": 10,
                "gwp": float(levels.loc[target_year]),
                "percentile": percentile,
                "mc_standard_error": np.sqrt(
                    percentile * (1.0 - percentile) / 10000
                ),
            }
        )
        q_start = fitted.x
    return pd.DataFrame(rows)


def calculate_annual_one_step(q_start: np.ndarray) -> pd.DataFrame:
    """Calculate annual one-step forecasts with yearly refitting."""
    updated = model_frame(include_update=True)
    rows = []
    for target_year in range(2020, 2026):
        prior = updated.loc[updated["year"] < target_year].copy()
        observed = updated.loc[updated["year"] == target_year].iloc[0]
        fitted = fit(prior, q_start)
        percentile = simulate_percentile(
            fitted,
            prior,
            start_value=float(prior.iloc[-1]["GWP"]),
            observed_value=float(observed["GWP"]),
            horizon=float(target_year - prior.iloc[-1]["year"]),
            seed=902381476 + target_year,
        )
        rows.append(
            {
                "year": target_year,
                "fit_through": target_year - 1,
                "horizon": 1,
                "gwp": float(observed["GWP"]),
                "percentile": percentile,
                "mc_standard_error": np.sqrt(
                    percentile * (1.0 - percentile) / 10000
                ),
            }
        )
        q_start = fitted.x
    return pd.DataFrame(rows)


def calculate_fixed_2010() -> pd.DataFrame:
    """Calculate forecasts that preserve the 2010 information set."""
    historical = model_frame(include_update=False)
    updated = model_frame(include_update=True)
    prior = historical.loc[historical["year"] <= 2010].copy()
    fitted = fit(prior)
    start_value = float(prior.iloc[-1]["GWP"])

    rows = []
    for target_year in range(2020, 2026):
        observed = updated.loc[updated["year"] == target_year].iloc[0]
        percentile = simulate_percentile(
            fitted,
            prior,
            start_value=start_value,
            observed_value=float(observed["GWP"]),
            horizon=float(target_year - 2010),
            seed=902381476,
        )
        rows.append(
            {
                "year": target_year,
                "fit_through": 2010,
                "horizon": target_year - 2010,
                "gwp": float(observed["GWP"]),
                "percentile": percentile,
                "mc_standard_error": np.sqrt(
                    percentile * (1.0 - percentile) / 10000
                ),
            }
        )
    return pd.DataFrame(rows)


def render_historical_svg(update: pd.DataFrame, output: Path) -> None:
    width, height = 1440, 1047
    left, right, top, bottom = 12, 125, 88, 85
    plot_width = width - left - right
    plot_height = height - top - bottom
    teal, red, ink = "#309ba6", "#cf3d32", "#2e2d2c"
    x_max, x_min = np.log(2035 - 1600), 1.75

    def sx(year: float) -> float:
        value = np.log(2035.0 - year)
        return left + (x_max - value) / (x_max - x_min) * plot_width

    def sy(percentile: float) -> float:
        return top + (1.0 - percentile) * plot_height

    def path(rows) -> str:
        return " ".join(
            ("M" if index == 0 else "L")
            + f" {sx(float(year)):.2f} {sy(float(percentile)):.2f}"
            for index, (year, percentile) in enumerate(rows)
        )

    grid = []
    for index in range(101):
        value = index / 100
        major = index % 10 == 0
        grid.append(
            f'<line x1="{left}" y1="{sy(value):.2f}" '
            f'x2="{left + plot_width}" y2="{sy(value):.2f}" '
            f'class="{"grid-major" if major else "grid-minor"}"/>'
        )
    ticks = "".join(
        f'<text x="{left + plot_width + 15}" y="{sy(value) + 9:.2f}" '
        f'class="tick">{value * 100:.0f}%</text>'
        for value in np.arange(0.0, 1.01, 0.1)
    )

    historical = list(PUBLISHED_PERCENTILES.items())
    historical_path = path(historical)
    update_path = path(
        [(2019, PUBLISHED_PERCENTILES[2019])]
        + list(update[["year", "percentile"]].itertuples(index=False, name=None))
    )
    history_points = "".join(
        f'<circle cx="{sx(year):.2f}" cy="{sy(value):.2f}" r="4.8" '
        f'class="historical-point"/>'
        for year, value in historical
    )
    update_points = "".join(
        f'<circle cx="{sx(row.year):.2f}" cy="{sy(row.percentile):.2f}" '
        f'r="5.4" class="update-point"/>'
        for row in update.itertuples(index=False)
    )

    historical_labels = "".join(
        f'<text x="{sx(year) + 7:.2f}" y="{sy(value) - 5:.2f}" '
        f'class="historical-label">{year}</text>'
        for year, value in historical
    )
    offsets = {
        2020: (8, 27, "start"),
        2021: (7, -13, "start"),
        2022: (-8, -19, "end"),
        2023: (8, 33, "start"),
        2024: (10, -22, "start"),
        2025: (9, 18, "start"),
    }
    update_labels = "".join(
        f'<text x="{sx(row.year) + offsets[row.year][0]:.2f}" '
        f'y="{sy(row.percentile) + offsets[row.year][1]:.2f}" '
        f'text-anchor="{offsets[row.year][2]}" class="update-label">'
        f'{row.year}</text>'
        for row in update.itertuples(index=False)
    )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
  text {{ font-family: Georgia, "Times New Roman", serif; fill: {ink}; }}
  .title {{ font-size: 36px; }}
  .axis-label {{ font-size: 34px; }}
  .tick {{ font-size: 26px; }}
  .historical-label {{ font-size: 24px; }}
  .update-label {{ font-size: 20px; fill: {red}; }}
  .grid-minor {{ stroke: #eeeeee; stroke-width: 1; }}
  .grid-major {{ stroke: #d7d7d7; stroke-width: 1.25; }}
  .historical-line {{ fill: none; stroke: {teal}; stroke-width: 2.5; }}
  .update-line {{ fill: none; stroke: {red}; stroke-width: 2.5; }}
  .historical-point {{ fill: {teal}; }}
  .update-point {{ fill: {red}; stroke: white; stroke-width: 1; }}
</style>
<rect width="{width}" height="{height}" fill="white"/>
<text x="{left}" y="49" class="title">{escape("Percentile of GWP in distribution when model fit to previous data")}</text>
<g>{"".join(grid)}</g>
<line x1="{left + plot_width}" y1="{top}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="{ink}" stroke-width="1.5"/>
<path d="{historical_path}" class="historical-line"/>
<path d="{update_path}" class="update-line"/>
{history_points}
{update_points}
{historical_labels}
{update_labels}
{ticks}
<text x="{left + plot_width / 2:.2f}" y="{height - 24}" text-anchor="middle" class="axis-label">Year</text>
</svg>
"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg)


def main() -> None:
    historical = model_frame(include_update=False)
    full_fit, standard_errors, percentile_2019 = validate_port(historical)
    rolling = calculate_rolling_decade(full_fit.x.copy())
    annual = calculate_annual_one_step(full_fit.x.copy())
    fixed = calculate_fixed_2010()
    render_historical_svg(rolling, MAIN_OUTPUT)
    render_historical_svg(annual, ANNUAL_ALTERNATIVE_OUTPUT)
    render_historical_svg(fixed, FIXED_ALTERNATIVE_OUTPUT)

    print("Validated full-sample estimates:")
    print("  theta:", np.array2string(unpack(full_fit.x), precision=8))
    print("  standard errors:", np.array2string(standard_errors, precision=8))
    print(f"Validated 2019 holdout percentile: {percentile_2019:.2%}")
    formatter = {
        "gwp": "{:,.3f}".format,
        "percentile": "{:.2%}".format,
        "mc_standard_error": "{:.2%}".format,
    }
    for label, result in (
        ("Rolling ten-year forecasts", rolling),
        ("Annual one-step alternative", annual),
        ("Fixed-2010 alternative", fixed),
    ):
        print(f"\n{label}:")
        print(result.to_string(index=False, formatters=formatter))
    print(f"\nWrote {MAIN_OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {ANNUAL_ALTERNATIVE_OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {FIXED_ALTERNATIVE_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
