## Summary

This draft extends the article's **"Percentile of GWP in distribution when
model fit to previous data"** chart from 2019 through 2025.

The new script ports the univariate Feller-diffusion likelihood and simulation
used by `Model GWP.do`/`asdf`, validates that port against archived results, and
then calculates rolling ten-year forecast percentiles for six new annual
observations.

## Before and after

| Published chart | Updated draft |
| --- | --- |
| ![Published prediction-percentile chart](https://coefficientgiving.org/wp-content/uploads/BernouDiffPredGWP12KDecBlog.png) | ![Prediction-percentile chart through 2025](https://raw.githubusercontent.com/llj0824/Modeling-Human-Trajectory/agent/update-gwp-percentiles-2025/reference-output-sample/gwp-prediction-percentiles-through-2025.svg) |

The published historical dots remain fixed at their displayed values. Red
identifies the added 2020-2025 rolling ten-year observations.

## How to read each point

Each historical dot is one rolling out-of-sample check:

1. Fit the model to all GWP observations through the preceding plotted year.
2. Start every simulated path at that preceding year's actual GWP.
3. Simulate 10,000 paths to the target year, drawing model parameters from
   their estimated covariance matrix.
4. Rank the target year's actual GWP among the 10,000 simulated endpoints.

The historical sequence means:

| Dot | Fit data through | Simulate |
| --- | --- | --- |
| 1600 | 1500 | 1500-1600 |
| 1700 | 1600 | 1600-1700 |
| 1820 | 1700 | 1700-1820 |
| 1870 | 1820 | 1820-1870 |
| 1913 | 1870 | 1870-1913 |
| 1940 | 1913 | 1913-1940 |
| 1950 | 1940 | 1940-1950 |
| 1960 | 1950 | 1950-1960 |
| 1970 | 1960 | 1960-1970 |
| 1980 | 1970 | 1970-1980 |
| 1990 | 1980 | 1980-1990 |
| 2000 | 1990 | 1990-2000 |
| 2010 | 2000 | 2000-2010 |
| 2019 | 2010 | 2010-2019 |

The added dots hold the forecast horizon constant at ten years:

| Dot | Fit data through | Simulate | Actual GWP percentile |
| --- | --- | --- | ---: |
| 2020 | 2010 | 2010-2020 | 13.87% |
| 2021 | 2011 | 2011-2021 | 15.42% |
| 2022 | 2012 | 2012-2022 | 15.32% |
| 2023 | 2013 | 2013-2023 | 15.79% |
| 2024 | 2014 | 2014-2024 | 15.35% |
| 2025 | 2015 | 2015-2025 | 14.61% |

Adjacent windows share nine of their ten years. The six dots should therefore
be read as several views of one sustained recent growth episode, not as six
independent observations.

A 15% dot means actual GWP was greater than about 15% of the simulated
endpoints and lower than about 85%. It does not mean the model assigned that
year a 15% probability.

## Data and method

The port reproduces the archived parameter estimates:

| Parameter | Archived Stata estimate | Python reproduction |
| --- | ---: | ---: |
| `ln(a)` | -12.66 | -12.661317 |
| `b` | 0.0000186 | 0.000018571 |
| `nu` | -23.78 | -23.7788 |
| `gamma` | -1.813 | -1.81289 |

- Historical GWP and uncertainty weights are reconstructed from `GWP.xlsx`
  using the same rules as the first `PrepData` pass in `Model GWP.do`.
- The update uses World Bank indicator `NY.GDP.MKTP.PP.KD`, GDP at constant
  2021 PPP, retrieved July 27, 2026.
- Each forecast uses 10,000 paths and 10,000 Euler steps, including parameter
  uncertainty.

## Validation

Before generating the new dots, the script reproduces the 2019 rolling
forecast at **20.71%**, consistent with the published dot at approximately 21%.

Monte Carlo standard errors for the six main percentiles are 0.35-0.36
percentage points.

## Alternative forecast constructions

The main chart uses a constant ten-year horizon. Two alternatives answer
different questions:

| Fixed 2010 information set | Annual one-step refits |
| --- | --- |
| ![Fixed-2010 prediction percentiles](https://raw.githubusercontent.com/llj0824/Modeling-Human-Trajectory/agent/update-gwp-percentiles-2025/reference-output-sample/gwp-prediction-percentiles-alternative-fixed-2010.svg) | ![Annual one-step prediction percentiles](https://raw.githubusercontent.com/llj0824/Modeling-Human-Trajectory/agent/update-gwp-percentiles-2025/reference-output-sample/gwp-prediction-percentiles-alternative-annual.svg) |

### Fixed 2010 information set

This version asks how each eventual outcome compares with the distribution
implied by information available in 2010. The model is fitted once, every path
starts at actual 2010 GWP, and the simulation horizon lengthens with each
target.

| Dot | Fit data through | Simulate | Actual GWP percentile |
| --- | --- | --- | ---: |
| 2020 | 2010 | 2010-2020 | 14.97% |
| 2021 | 2010 | 2010-2021 | 15.75% |
| 2022 | 2010 | 2010-2022 | 14.57% |
| 2023 | 2010 | 2010-2023 | 13.53% |
| 2024 | 2010 | 2010-2024 | 12.66% |
| 2025 | 2010 | 2010-2025 | 11.71% |

### Annual one-step refits

This version asks how the next year's outcome compares with a model refitted
through the immediately preceding year.

| Dot | Fit data through | Simulate | Actual GWP percentile |
| --- | --- | --- | ---: |
| 2020 | 2019 | 2019-2020 | 18.07% |
| 2021 | 2020 | 2020-2021 | 50.54% |
| 2022 | 2021 | 2021-2022 | 39.57% |
| 2023 | 2022 | 2022-2023 | 38.16% |
| 2024 | 2023 | 2023-2024 | 37.91% |
| 2025 | 2024 | 2024-2025 | 38.44% |

The alternatives are included for comparison and are not used for the red dots
in the main updated chart.

## Deliberate limitations

- Annual starting observations from 2011-2015 make the main extension a new
  rolling-horizon diagnostic rather than a literal continuation of the paper's
  preferred decennial sample.
- The overlapping ten-year windows are not statistically independent.
- The percentiles locate actual GWP within model-generated distributions; they
  do not identify the causes of slower growth.
- This PR does not revise the paper's median takeoff year, probability of no
  takeoff, or other stochastic-model estimates.
