# Optimise smoothing parameter for a given tree and branch statistic

This function maximises the deviance explained of future evolution of
the given branch statistic via its dependence on a smoothing parameter
(tau). This function must be given a tree and a function (e.g. \`codls\`
or \`lbi\`) which takes a tree and a smoothing parameter and returns a
given branch statistic. The association of the branch statistic and
future coalescent events (prediction of evolution) is modelled using a
GAM and adjusting for branch lengths and distance from the root.

## Usage

``` r
optsmooth(
  tr,
  func,
  logtaulb = -4,
  logtauub = 35,
  startpc = 50,
  endpc = 100,
  nobj = 100
)
```

## Arguments

- tr:

  A phylogenetic tree in ape::phylo format

- logtaulb:

  Lower bound of precision parameteters

- logtauub:

  Upper bound of precision parameteters

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

## Value

Output of \`optimise\`. \$minimum contains the optimal smoothing
parameter
