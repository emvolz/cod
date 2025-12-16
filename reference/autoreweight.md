# Automatically reweight sample units that may be oversampled

Given a fit from \`codls\` and a set of samples which are suspected of
over-sampling, this function will re-compute \`codls\` over a range of
reweighted samples. This will identify the sample weight at which an
association is lost between coalescent odds (psi) and the given set of
samples. This is an appropriate weight to use if there is an association
between coalescent odds and psi that is due to sampling effects and not
due to evolutionary effects, but note that this method may mask
evolutionary effects if any are present.

## Usage

``` r
autoreweight(f, rwtips, wlb = 0.01, wub = 0.5, res = 10, alpha = 0.05)
```

## Arguments

- f:

  A \`codls\` fit

- rwtips:

  Vector of samples (type character) which are suspected of
  over-sampling

- wlb:

  Numeric lower bound of sample weights to examine

- wub:

  Numeric upper bound of sample weights to examine

- res:

  Integer number of weights to examine

- alpha:

  The p value threshold used for selecting the optimal weight

## Value

A list with components \`fit\`: the reweigthed \`codls\` fit;
\`weights\`: a new vector of sample weights; \`optimalweight\` the
scalar weight applied to oversampled units; and \`summary\`: a data
frame showing regression p values over a range of weights
