# Fit a genealogical placement GMRF model using maximum likelihood

Fits the COD GMRF model using the \`mgcv::gam\` method. Additional
arguments can be passed to \`gam\`; see documentation for that method.
Using method="REML" can speed execution by using a constrained maximum
likelihood approach. Additionally, an approximate reduced-rank MRF model
can be fitted by supplying the \`k\` parameter. If tau is not provided,
\`codls\` is also used to optimise this parameter. This method is slower
than \`codls\` and is not recommended for trees with more than several
hundred samples.

## Usage

``` r
codbinomial(tr1, logtau = NULL, k = Inf, profcontrol = list(), ...)
```

## Arguments

- tr1:

  Phylogenetic tree in ape::phylo format

- logtau:

  Precision parameter. If NULL, will invoke \`tauprofile\` to find best
  value.

- k:

  Fits a reduced-rank MRF model if k is an integer \< number of nodes in
  the input tree. This can speed calculation but reduces precision.

- profcontrol:

  Optional list of arguments passed to \`tauprofile\`

- ...:

  Additional arguments are passed to \`mgcv::gam\`

## Value

A COD GMRF model fit. Includes the GAM model fit.

## Details

The ML COD GMRF method does not currently support inverse probability
weighting of samples. Use \`codls\` if sample weighting is needed.
