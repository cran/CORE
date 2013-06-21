\name{PCOREiteration}
\alias{PCOREiteration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PCOREiteration(z, maxmark, minscore = 0, nprof = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{maxmark}{
%%     ~~Describe \code{maxmark} here~~
}
  \item{minscore}{
%%     ~~Describe \code{minscore} here~~
}
  \item{nprof}{
%%     ~~Describe \code{nprof} here~~
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
function (z, maxmark, minscore = 0, nprof = 1) 
{
    chu <- unique(z[, "chrom"])
    chc <- rep(0, length(chu))
    for (i in 1:length(chu)) chc[i] <- sum(z[, "chrom"] == chu[i])
    for (ich in 1:length(chu)) {
        za <- z[z[, "chrom"] == chu[rev(order(chc))][ich], c("start", 
            "end", "weight"), drop = F]
        za <- za[order(za[, "start"]), , drop = F]
        winloci <- matrix(ncol = 4, nrow = maxmark, dimnames = list(NULL, 
            c("chrom", "start", "end", "score")))
        winloci[, "chrom"] <- chu[rev(order(chc))][ich]
        for (i in 1:maxmark) {
            if (nrow(za) == 0) 
                break
            y <- cbind(c(za[, "end"] + 1, za[, "start"]), c(-za[, 
                "weight"], za[, "weight"]))
            y <- y[order(y[, 1]), , drop = F]
            cy2 <- cumsum(y[, 2])
            score <- max(cy2)
            if (score < minscore) 
                break
            stabstart <- y[which.max(cy2), 1]
            stabend <- min(y[y[, 1] > stabstart & y[, 2] <= 0, 
                1]) - 1
            winloci[i, c("start", "end", "score")] <- c(stabstart, 
                stabend, score)
            if (ich > 1) 
                if (((sum(accu[, "score"] > winloci[i, "score"]) + 
                  i) >= maxmark)) 
                  break
            za <- za[!(za[, "start"] <= stabstart & za[, "end"] >= 
                stabstart), , drop = F]
        }
        winloci <- winloci[!is.na(winloci[, "score"]), , drop = F]
        if (ich == 1) 
            accu <- winloci
        if (ich > 1) {
            accu <- rbind(accu, winloci)
            accu <- accu[accu[, "score"] >= minscore, , drop = F]
            accu <- accu[order(accu[, "score"], decreasing = T), 
                , drop = F][1:min(maxmark, nrow(accu)), , drop = F]
        }
    }
    accu[, "score"] <- accu[, "score"]/nprof
    return(accu)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
