\name{fitgpgmrf}
\alias{fitgpgmrf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitgpgmrf(tr1, tau = c(1, NULL), profcontrol = list(), inverseprobabilityweights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr1}{
%%     ~~Describe \code{tr1} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{profcontrol}{
%%     ~~Describe \code{profcontrol} here~~
}
  \item{inverseprobabilityweights}{
%%     ~~Describe \code{inverseprobabilityweights} here~~
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
function (tr1, tau = c(1, NULL), profcontrol = list(), inverseprobabilityweights = NULL) 
{
    if (!inherits(tr1, "cggephylo") & inherits(tr1, "phylo")) {
        tr1 <- .maketreedata(tr1)
    }
    tau = tau[1]
    tpdf <- NULL
    if (is.null(tau)) {
        tpargs <- modifyList(TPARGS, profcontrol)
        tpargs$tr = tr1
        tpargs$ipw = inverseprobabilityweights
        tpdf <- do.call(tauprofile, tpargs)
        tau <- tpdf$tau[which.min(tpdf$loss)]
        print(tpdf)
    }
    st1 <- Sys.time()
    whno = tr1$whno
    nr = tr1$nr
    i <- which(!is.na(tr1$parent))
    ary = rep(0, nr)
    y = c(ary, tr1$nodey)
    arw = tau/tr1$brlens
    nodew <- .computenodew(tr1, inverseprobabilityweights)
    w = c(arw, nodew)
    st2b <- Sys.time()
    ai = 1:nr
    aj = whno
    ax = rep(1, nr)
    ai = c(ai, 1:nr)
    aj = c(aj, tr1$parent[whno])
    ax = c(ax, rep(-1, nr))
    k <- nr + 1
    coindices <- list()
    coii <- 1
    for (co in tr1$coalescentcohorts) {
        coindices[[coii]] <- k:(k + length(co) - 1)
        coii <- coii + 1
        k <- k + length(co)
    }
    ncoi <- sum(sapply(tr1$coalescentcohorts, length))
    coi <- (nr + 1):(nr + ncoi)
    coj <- do.call(c, tr1$coalescentcohorts)
    cox <- rep(1, ncoi)
    ai <- c(ai, coi)
    aj <- c(aj, coj)
    ax <- c(ax, cox)
    X <- Matrix::sparseMatrix(i = ai, j = aj, x = ax)
    W <- Matrix::Diagonal(x = w)
    QQ <- t(X) \%*\% W \%*\% X
    b <- t(X) \%*\% W \%*\% y
    f2beta = as.vector(solve(QQ, b))
    st3 <- Sys.time()
    structure(list(coef = f2beta, tau = tau, data = tr1, X = X, 
        W = W, y = y, nr = nr, arindices = 1:nr, logoddsindices = (nr + 
            1):nrow(X), istartnodeterms = nr + 1, coindices = coindices, 
        tauprofile = tpdf, runtime = st3 - st1), class = "gpgmrf")
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
