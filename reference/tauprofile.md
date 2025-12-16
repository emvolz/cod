# Evaluate the loss function of the cod model across a range of tau (precision parameter) values

Evaluate the loss function of the cod model across a range of tau
(precision parameter) values

## Usage

``` r
tauprofile(
  tr,
  logtaulb = -4,
  logtauub = 35,
  res = 11,
  startpc = 75,
  endpc = 100,
  nobj = 100,
  ipw = NULL,
  ncpu = 1
)
```

## Arguments

- tr:

  A phylogenetic tree in ape::phylo format

- logtaulb:

  Lower bound of precision parameteters

- logtauub:

  Upper bound of precision parameteters

- res:

  Number of tau values to evaluate

- startpc:

  The initial per cent of nodes in the tree counting from root to tips
  where the loss function will be evaluated

- endpc:

  The final per cent of nodes in the tree counting from root to tips
  where the loss function will be evaluated

- nobj:

  The integer number of points along the tree where the loss function
  will be evaluated. If Inf, will use all points between startpc and
  endpc, but may be slow.

- ipw:

  Optional inverse probability weights for each sample

## Value

A data frame containing the loss function evaluated over a range of tau
values
