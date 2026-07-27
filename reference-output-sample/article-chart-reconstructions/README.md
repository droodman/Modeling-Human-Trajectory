# Reconstructed article charts

These three SVGs rebuild the article's introductory visual sequence from the
author's `GWP.xlsx` data and the core `PrepData` transformation in `Model GWP.do`:

1. ordinary axes;
2. a logarithmic vertical axis plus an exponential fit; and
3. logarithmic axes after transforming time to "years until 2047," plus a
   power-law fit.

They are independent visual reconstructions, not image copies and not a
rerun of the paper's full statistical model. The source code is
`scripts/reconstruct_article_charts.py`.
# Article-chart reconstructions

`01-ordinary-axes.svg`, `02-log-y-exponential-fit.svg`, and
`03-transformed-time-power-fit.svg` reconstruct the paper-era series ending in
2019. Their marks are blue.

The matching `*-through-2025.svg` files add a visually distinct red layer for
2020–2025. Those points use World Bank world PPP GDP growth, chain-linked to
the reconstructed 2019 value. They are a chart extension, not a re-estimation
of the paper's stochastic model.
