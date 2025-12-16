# A Capitalized Title (ideally limited to 65 characters)

## Usage

``` r
coef.gpgmrf(f)
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
    f$coef
  }
#> function (f) 
#> {
#>     stopifnot(inherits(f, "gpgmrf"))
#>     f$coef
#> }
#> <environment: 0x55d88410cd68>
```
