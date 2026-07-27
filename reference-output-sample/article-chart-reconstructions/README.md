# Article-chart comparisons: 2019 source and 2025 visual extension

## Review pairs

The **before** images are David Roodman's published 2019 charts. They remain
hosted at the article's original URLs; they are not files in this repository.
The **after** images below are programmatically generated local files in this
directory. Each keeps the reconstructed 2019 history in the article's teal and adds a red
2020–2025 segment.

| View | Before: published 2019 chart | After: generated 2025 extension |
| --- | --- | --- |
| Ordinary axes | [Open original PNG](https://coefficientgiving.org/wp-content/uploads/Roodman_GWP_10000_BCE-2019_1.png) | [Open SVG](01-ordinary-axes-through-2025.svg) · [Open PNG](01-ordinary-axes-through-2025.png) |
| Log y-axis / exponential-fit view | [Open original PNG](https://coefficientgiving.org/wp-content/uploads/Roodman_GWP_10000_BCE-2019_2.png) | [Open SVG](02-log-y-exponential-fit-through-2025.svg) · [Open PNG](02-log-y-exponential-fit-through-2025.png) |
| Transformed-time / power-fit view | [Open original PNG](https://coefficientgiving.org/wp-content/uploads/Roodman_GWP_10000_BCE-2019_3.png) | [Open SVG](03-transformed-time-power-fit-through-2025.svg) · [Open PNG](03-transformed-time-power-fit-through-2025.png) |

The local `01-ordinary-axes.svg`, `02-log-y-exponential-fit.svg`, and
`03-transformed-time-power-fit.svg` are generated 2019 reconstructions used
as the teal baseline. They are not copies of the published PNGs.

## How the after charts are made

These three SVGs rebuild the article's introductory visual sequence from the
author's `GWP.xlsx` data and the core `PrepData` transformation in `Model GWP.do`:

1. ordinary axes;
2. a logarithmic vertical axis plus an exponential fit; and
3. logarithmic axes after transforming time to "years until 2047," plus a
   power-law fit.

They are independent visual reconstructions, not image copies and not a
rerun of the paper's full statistical model. The source code is
`scripts/reconstruct_article_charts.py`.
## Color-coded extensions

`01-ordinary-axes.svg`, `02-log-y-exponential-fit.svg`, and
`03-transformed-time-power-fit.svg` reconstruct the paper-era series ending in
2019. Their marks use the article's teal.

The matching `*-through-2025.svg` files add a visually distinct red layer for
2020–2025. Those points use World Bank world PPP GDP growth, chain-linked to
the reconstructed 2019 value. They are a chart extension, not a re-estimation
of the paper's stochastic model.
