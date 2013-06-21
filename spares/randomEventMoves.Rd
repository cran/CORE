\name{randomEventMoves}
\alias{randomEventMoves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
randomEventMoves(eventLengths, boundaries)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eventLengths}{
%%     ~~Describe \code{eventLengths} here~~
}
  \item{boundaries}{
%%     ~~Describe \code{boundaries} here~~
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
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (eventLengths, boundaries) 
{
    chrlengths <- boundaries[, "end"] - boundaries[, "start"] + 
        1
    ochl <- order(chrlengths, decreasing = T)
    schl <- chrlengths[ochl]
    cschl <- cumsum(schl)
    sb <- boundaries[ochl, ]
    oe <- order(eventLengths, decreasing = T)
    z <- cbind(c(chrlengths, eventLengths), c(rep(1, length(chrlengths)), 
        rep(0, length(eventLengths))))
    z <- z[order(z[, 1], decreasing = T), ]
    chrmax <- cumsum(z[, 2])[z[, 2] == 0]
    eventLengths <- z[z[, 2] == 0, 1]
    abslen <- floor(cschl[chrmax] * runif(length(chrmax)))
    z <- cbind(c(cschl, abslen), c(rep(1, length(cschl)), rep(0, 
        length(abslen))), c(rep(0, length(cschl)), eventLengths))
    oz <- order(z[, 1])
    oze <- order(z[z[, 3] != 0, 1])
    z <- z[oz, ]
    z[, 2] <- cumsum(z[, 2]) + 1
    z <- z[z[, 3] != 0, , drop = F]
    start <- sb[z[, 2], "start"] + floor((schl[z[, 2]] - z[, 
        3] + 1) * runif(nrow(z)))
    end <- start + z[, 3] - 1
    chrom <- ochl[z[, 2]]
    sec <- matrix(ncol = 3, nrow = length(start), dimnames = list(NULL, 
        c("start", "end", "chrom")))
    sec[oze, ] <- cbind(start, end, chrom)
    sec[oe, ] <- sec
    return(sec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
