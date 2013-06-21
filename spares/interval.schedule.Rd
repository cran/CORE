\name{interval.schedule}
\alias{interval.schedule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interval.schedule(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
%%     ~~Describe \code{v} here~~
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
function (v) 
{
    if (!("start" \%in\% dimnames(v)[[2]]) | !("end" \%in\% dimnames(v)[[2]]) | 
        length(dim(v)) != 2) 
        stop("interval.schedule: incorrect form of input\n")
    lv <- nrow(v)
    if (sum(sort(v[, "start"]) == v[, "start"]) != lv) 
        stop("interval.schedule: input incorrectly ordered\n")
    unsch <- lv
    schedule <- 0
    vsched <- rep(0, lv)
    while (unsch > 0) {
        schedule <- schedule + 1
        vsched[cummax(v[, "end"] * (vsched == 0)) == v[, "end"]] <- schedule
        unsch <- sum(vsched == 0)
    }
    return(vsched)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
