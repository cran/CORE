\name{fill.values}
\alias{fill.values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fill.values(fromto, values, ldest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fromto}{
%%     ~~Describe \code{fromto} here~~
}
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{ldest}{
%%     ~~Describe \code{ldest} here~~
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
function (fromto, values, ldest) 
{
    dest <- rep(0, ldest)
    z <- cbind(c(fromto[, 1], fromto[, 2] + 1), c(values, -values))
    z <- z[order(z[, 1]), ]
    z[, 2] <- cumsum(z[, 2])
    z <- z[nrow(z) - rev(match(unique(rev(z[, 1])), rev(z[, 1]))) + 
        1, ]
    if (z[nrow(z), 1] > ldest) 
        z <- matrix(ncol = ncol(z), data = z[-nrow(z), ])
    dest[z[, 1]] <- z[, 2]
    zz <- rep(0, length(dest))
    zz[z[, 1]] <- c(0, z[-nrow(z), 2])
    return(cumsum(dest - zz))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
