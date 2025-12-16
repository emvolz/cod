# Root mean square coalescent log odds

This is a summary statistic that describes the amount of variation in
coalescent rates across lineages in a phylogenetic tree It is defined as
\$\$ \sqrt{ \sum_i l_i \psi_i^2 / L } \$\$ where the sum is over all
branches i in the tree and weighted by branch length \\l_i\\ and where
\\ L = \sum_i l_i \\

## Usage

``` r
rmsclo(f)
```

## Arguments

- f:

  Fit from \`codls\`

## Value

Numeric RMSCLO
