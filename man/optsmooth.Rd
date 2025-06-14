% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cppe2.R
\name{optsmooth}
\alias{optsmooth}
\title{Optimise smoothing parameter for a given tree and branch statistic}
\usage{
optsmooth(
  tr,
  func,
  logtaulb = -4,
  logtauub = 35,
  startpc = 50,
  endpc = 100,
  nobj = 100
)
}
\arguments{
\item{tr}{A phylogenetic tree in ape::phylo format}

\item{logtaulb}{Lower bound of precision parameteters}

\item{logtauub}{Upper bound of precision parameteters}

\item{startpc}{The initial per cent of nodes in the tree counting from root to tips where the loss function will be evaluated}

\item{endpc}{The final per cent of nodes in the tree counting from root to tips where the loss function will be evaluated}

\item{nobj}{The integer number of points along the tree where the loss function will be evaluated. If Inf, will use all points between startpc and endpc, but may be slow.}
}
\value{
Output of `optimise`. $minimum contains the optimal smoothing parameter
}
\description{
The minimum minimum prediction deviance is the criterian being maximised. 
This function must be given a tree and a function (e.g. `codls` or `lbi`) which takes a tree and a smoothing parameter and returns a given branch statistic. 
The branch statistic is associated with coalescence with future lineages (prediction of evolution) using a GAM and adjusting for branch lengths and distance from the root.
}
