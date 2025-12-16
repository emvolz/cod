# Compute CH index over a range of clustering thresholds

Compute CH index over a range of clustering thresholds

## Usage

``` r
chindices(f, clths = seq(0.1, 1.5, length = 20), rescale = TRUE)
```

## Arguments

- f:

  A fit from \`codls\`

- clths:

  Numeric vector of clustering thresholds

- rescale:

  Passed to \`computeclusters\`
