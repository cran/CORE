\name{makepairs}
\alias{makepairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makepairs(events)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{events}{
%%     ~~Describe \code{events} here~~
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
function (events) 
{
    if ("chrom" \%in\% dimnames(events)[[2]]) {
        events <- events[order(events[, "chrom"]), , drop = F]
        chromstarts <- tapply(X = events[, "start"], INDEX = events[, 
            "chrom"], FUN = min)
        chromends <- tapply(X = events[, "end"], INDEX = events[, 
            "chrom"], FUN = max)
    }
    else {
        chromstarts <- min(events[, "start"])
        chromends <- max(events[, "end"])
    }
    ustarts <- sort(unique(events[, "start"]))
    uends <- sort(unique(events[, "end"]))
    z <- cbind(c(chromstarts, ustarts), c(rep(1, length(chromstarts)), 
        rep(0, length(ustarts))), c(rep(0, length(chromstarts)), 
        1:length(ustarts)))
    z <- z[order(z[, 1]), ]
    sch <- cumsum(z[, 2])[z[, 3] != 0][order(z[z[, 3] != 0, 3])]
    z <- cbind(c(chromstarts, uends), c(rep(1, length(chromstarts)), 
        rep(0, length(uends))), c(rep(0, length(chromstarts)), 
        1:length(uends)))
    z <- z[order(z[, 1]), ]
    ech <- cumsum(z[, 2])[z[, 3] != 0][order(z[z[, 3] != 0, 3])]
    sfirst <- match(1:length(chromstarts), sch)[ech]
    z <- cbind(c(ustarts, uends), c(rep(1, length(ustarts)), 
        rep(0, length(uends))))
    z <- z[order(z[, 1]), ]
    slast <- cumsum(z[, 2])[z[, 2] == 0]
    z <- cbind(rep(1:length(uends), times = (slast - sfirst + 
        1)), rep(0, sum(slast - sfirst + 1)))
    z[match(1:length(uends), z[, 1]), 2] <- match(1:length(uends), 
        z[, 1])
    mypairs <- matrix(ncol = 3, dimnames = list(NULL, c("start", 
        "end", "id")), data = c(ustarts[(1:nrow(z)) - cummax(z[, 
        2]) + sfirst[z[, 1]]], uends[z[, 1]], 1:length(uends[z[, 
        1]])))
    return(mypairs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
