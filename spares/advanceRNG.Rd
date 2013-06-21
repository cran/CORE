\name{advanceRNG}
\alias{advanceRNG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
advanceRNG(randopt = c("SIMPLE", "RESCALE"), nrand, nevents)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{randopt}{
%%     ~~Describe \code{randopt} here~~
}
  \item{nrand}{
%%     ~~Describe \code{nrand} here~~
}
  \item{nevents}{
%%     ~~Describe \code{nevents} here~~
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
function (randopt = c("SIMPLE", "RESCALE"), nrand, nevents) 
{
    randopt <- match.arg(randopt)
    switch(randopt, SIMPLE = replicate(nrand, {
        a <- runif(2 * nevents)
        rm(a)
    }), RESCALE = replicate(nrand, {
        a <- sample(nevents)
        a <- runif(nevents)
        rm(a)
    }))
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
