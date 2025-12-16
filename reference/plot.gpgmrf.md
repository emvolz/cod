# A Capitalized Title (ideally limited to 65 characters)

## Usage

``` r
plot.gpgmrf(f)
```

## Arguments

- f:

## Details

## Value

## References

## Author

## Note

## See also

## Examples

``` r
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--  or standard data sets, see data().

## The function is currently defined as
function (f) 
{
    stopifnot(inherits(f, "gpgmrf"))
    f2beta = f$coef
    tr1 = f$data
    class(tr1) <- "phylo"
    fdf <- data.frame(node = 1:length(tr1$nodetimes), theta = f2beta)
    gtr1 = ggtree::ggtree(tr1) %<+% fdf
    gtr1 + aes(color = theta) + scale_color_gradient2(low = "blue", 
        mid = "lightblue", high = "red", midpoint = 0, limits = range(fdf$theta), 
        name = "ψ") + ggtree::theme_tree2() + ggtree::geom_tiplab()
  }
#> function (f) 
#> {
#>     stopifnot(inherits(f, "gpgmrf"))
#>     f2beta = f$coef
#>     tr1 = f$data
#>     class(tr1) <- "phylo"
#>     fdf <- data.frame(node = 1:length(tr1$nodetimes), theta = f2beta)
#>     gtr1 = ggtree::ggtree(tr1) %<+% fdf
#>     gtr1 + aes(color = theta) + scale_color_gradient2(low = "blue", 
#>         mid = "lightblue", high = "red", midpoint = 0, limits = range(fdf$theta), 
#>         name = "ψ") + ggtree::theme_tree2() + ggtree::geom_tiplab()
#> }
#> <environment: 0x55d880f922e8>
```
