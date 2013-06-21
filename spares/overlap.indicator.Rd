\name{overlap.indicator}
\alias{overlap.indicator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overlap.indicator(vstart, vend, wstart, wend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vstart}{
%%     ~~Describe \code{vstart} here~~
}
  \item{vend}{
%%     ~~Describe \code{vend} here~~
}
  \item{wstart}{
%%     ~~Describe \code{wstart} here~~
}
  \item{wend}{
%%     ~~Describe \code{wend} here~~
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
function (vstart, vend, wstart, wend) 
{
    lw <- length(wstart)
    lv <- length(vstart)
    z <- cbind(c(wstart, vend), c(rep(0, lw), 1:lv), c(1:lw, 
        rep(0, lv)))
    z <- z[order(z[, 1]), ]
    endbefore <- cummax(z[, 2])[order(z[, 3])][sort(z[, 3]) != 
        0]
    z <- cbind(c(vstart, wend), c(1:lv, rep((lv + 1), lw)), c(rep(0, 
        lv), 1:lw))
    z <- z[order(z[, 1]), ]
    startafter <- rev(cummin(rev(z[, 2])))[order(z[, 3])][sort(z[, 
        3]) != 0]
    return(cbind(endbefore + 1, startafter - 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
