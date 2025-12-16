# Compute phylogenetic clusters by cutting tree at branches with large changes in coalescent odds

Compute phylogenetic clusters by cutting tree at branches with large
changes in coalescent odds

## Usage

``` r
computeclusters(f, clth = NULL, rescale = TRUE, includeinternals = TRUE)
```

## Arguments

- f:

  A model fit from \`codls\`

- clth:

  Numeric threshold change in coalescent log odds. If NULL, will guess a
  good threshold based on the CH index

- rescale:

  if TRUE (default), coalescent log odds are rescaled (mean zero, unit
  variance) prior to applying thresholds

- includeinternals:

  if TRUE (default), internal nodes are also included in cluster output

## Value

A data frame with cluster asignment for each tip
