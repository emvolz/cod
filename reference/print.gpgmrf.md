# A Capitalized Title (ideally limited to 65 characters)

## Usage

``` r
print.gpgmrf(f)
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
    cat(" Genealogical placement GMRF model fit \n")
    print(f$data)
    cat("Range of coefficients: \n")
    print(range(coef(f)))
    cat(glue::glue("Precision parameter (tau): {f$tau} \n"))
    cat("\n")
    invisible(f)
  }
#> function (f) 
#> {
#>     stopifnot(inherits(f, "gpgmrf"))
#>     cat(" Genealogical placement GMRF model fit \n")
#>     print(f$data)
#>     cat("Range of coefficients: \n")
#>     print(range(coef(f)))
#>     cat(glue::glue("Precision parameter (tau): {f$tau} \n"))
#>     cat("\n")
#>     invisible(f)
#> }
#> <environment: 0x55d87ccefd70>
```
