# Fit a genealogical placement GMRF model using maximum likelihood

Fits the COD GMRF model using maximum likelihood. Additional arguments
are passed to \`optim\`. If tau is not provided, \`codls\` is also used
to optimise this parameter. This method is slower than \`codls\` and is
not recommended for trees with more than several hundred samples.

## Usage

``` r
codml(tr1, logtau = NULL, profcontrol = list(), ...)
```

## Arguments

- tr1:

  Phylogenetic tree in ape::phylo format

- logtau:

  Precision parameter. If NULL, will invoke \`tauprofile\` to find best
  value.

- profcontrol:

  Optional list of arguments passed to \`tauprofile\`

- ...:

  Additional arguments are passed to \`mgcv::gam\`

## Value

A COD GMRF model fit. Includes the GAM model fit.

## Details

The ML COD GMRF method does not currently support inverse probability
weighting of samples. Use \`codls\` if sample weighting is needed.
