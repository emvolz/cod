# Fit a COD GMRF model using weighted least squares

Fit a COD GMRF model using weighted least squares

## Usage

``` r
codls(tr1, logtau = NULL, profcontrol = list(), weights = NULL, ncpu = 1)
```

## Arguments

- tr1:

  Phylogenetic tree in ape::phylo format

- logtau:

  Precision parameter. If NULL, will invoke \`tauprofile\` to find best
  value.

- profcontrol:

  Optional list of arguments passed to \`tauprofile\`

- weights:

  An optional vector (named or unnamed) of sample weights for each tip
  in the input tree

- ncpu:

  Integer number of cpu's to use if parallel calculation of tau profile
  is desired

## Value

A COD GMRF model fit

## Examples

``` r
# A simple example that does not have population structure 
set.seed( 1111 )
tr <- ape::rtree( 100 )
f <- codls(tr)
#> c(-4, 0.0999999999999996, 4.2, 8.3, 12.4, 16.5, 20.6, 24.7, 28.8, 32.9, 37)c(839.610222117401, 372.398296997339, 345.306824047442, 345.280902945219, 345.28089580974, 345.28089580778, Inf, Inf, Inf, Inf, Inf)c("", "", "", "", "", "***", "", "", "", "", "")
#> Warning: Inf replaced by maximum positive value
#> Warning: Inf replaced by maximum positive value
#> Warning: Inf replaced by maximum positive value
#> c(-4, 0.0999999999999996, 4.2, 8.3, 12.4, 16.5, 20.6, 24.7, 28.8, 32.9, 37)c(839.610222117401, 372.398296997339, 345.306824047442, 345.280902945219, 345.28089580974, 345.28089580778, Inf, Inf, Inf, Inf, Inf)c("", "", "", "", "", "***", "", "", "", "", "")
coef(f) |> head() 
#> [1] -3.640017e-14 -3.132941e-14 -3.446806e-14 -3.171095e-14 -4.953729e-14
#> [6] -4.956475e-14
summary(f)
#>  Genealogical placement GMRF model fit 
#> 
#> Phylogenetic tree with 100 tips and 99 internal nodes.
#> 
#> Tip labels:
#>   t13, t36, t48, t79, t65, t66, ...
#> 
#> Rooted; includes branch length(s).
#> Range of coefficients: 
#> -7.19882172734318e-14 8.88092382524843e-14
#> Precision parameter (log tau): 16.761414287381 
#> 
#>   logprecision       RMSCLO Neff
#> 1     16.76141 4.240497e-14  100
#> 
if (FALSE) { # \dontrun{
plot(f)
} # }

# This example has population structure 
tr0 = rcoal(20); tr0$edge.length <- .01*tr0$edge.length 
tr1 = rcoal(80); 
dx <- (max(node.depth.edgelength( tr1 ))-max(node.depth.edgelength( tr0 )))
tr0$root.edge <- dx
tr <- bind.tree(tr0,tr1, position = dx)
f <- codls(tr)
#> c(-4, 0.0999999999999996, 4.2, 8.3, 12.4, 16.5, 20.6, 24.7, 28.8, 32.9, 37)c(3743.24780055663, 3210.78459384594, 3228.34397023542, 3228.61108342026, 3228.61307940198, Inf, 3228.569943925, 3228.569943925, Inf, Inf, Inf)c("", "***", "", "", "", "", "", "", "", "", "")
#> c(-4, 0.0999999999999996, 4.2, 8.3, 12.4, 16.5, 20.6, 24.7, 28.8, 32.9, 37)c(3743.24780055663, 3210.78459384594, 3228.34397023542, 3228.61108342026, 3228.61307940198, Inf, 3228.569943925, 3228.569943925, Inf, Inf, Inf)c("", "***", "", "", "", "", "", "", "", "", "")
summary(f) 
#>  Genealogical placement GMRF model fit 
#> 
#> Phylogenetic tree with 100 tips and 98 internal nodes.
#> 
#> Tip labels:
#>   t18, t1, t16, t9, t7, t11, ...
#> 
#> Rooted; includes branch length(s).
#> Range of coefficients: 
#> -0.548001884420184 1.15812937469376
#> Precision parameter (log tau): 0.83441648702799 
#> 
#>   logprecision    RMSCLO     Neff
#> 1    0.8344165 0.5697688 92.94934
#> 
if (FALSE) { # \dontrun{
plot(f)
} # }
```
