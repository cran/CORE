\name{randomRescaledEventMoves}
\alias{randomRescaledEventMoves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
randomRescaledEventMoves(events, boundaries, countOnly = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{events}{
%%     ~~Describe \code{events} here~~
}
  \item{boundaries}{
%%     ~~Describe \code{boundaries} here~~
}
  \item{countOnly}{
%%     ~~Describe \code{countOnly} here~~
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
function (events, boundaries, countOnly = F) 
{
    chrlengths <- boundaries[, "end"] - boundaries[, "start"] + 
        1
    chrfrac <- (events[, "end"] - events[, "start"] + 1)/chrlengths[match(events[, 
        "chrom"], boundaries[, "chrom"])]
    chrstart <- cumsum(chrlengths) - chrlengths + 1
    newchrom <- sample(unique(events[, "chrom"]), size = nrow(events), 
        replace = T)
    newstart <- chrstart[newchrom] + floor((1 - chrfrac) * runif(nrow(events)) * 
        chrlengths[match(newchrom, boundaries[, "chrom"])])
    if (countOnly) 
        return()
    newend <- newstart + pmax(0, round(chrfrac * chrlengths[match(newchrom, 
        boundaries[, "chrom"])]) - 1)
    return(matrix(ncol = 3, data = c(newstart, newend, newchrom), 
        dimnames = list(NULL, c("start", "end", "chrom"))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
