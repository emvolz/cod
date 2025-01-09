\name{plot.gpgmrf}
\alias{plot.gpgmrf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.gpgmrf(f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (f) 
{
    stopifnot(inherits(f, "gpgmrf"))
    f2beta = f$coef
    tr1 = f$data
    class(tr1) <- "phylo"
    fdf <- data.frame(node = 1:length(tr1$nodetimes), theta = f2beta)
    gtr1 = ggtree::ggtree(tr1) \%<+\% fdf
    gtr1 + aes(color = theta) + scale_color_gradient2(low = "blue", 
        mid = "lightblue", high = "red", midpoint = 0, limits = range(fdf$theta), 
        name = "ψ") + ggtree::theme_tree2() + ggtree::geom_tiplab()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
