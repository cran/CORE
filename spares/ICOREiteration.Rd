\name{ICOREiteration}
\alias{ICOREiteration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ICOREiteration(z, maxmark, pow = 1, minscore = 0, nprof = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{maxmark}{
%%     ~~Describe \code{maxmark} here~~
}
  \item{pow}{
%%     ~~Describe \code{pow} here~~
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
function (z, maxmark, pow = 1, minscore = 0, nprof = 1) 
{
    chu <- unique(z[, "chrom"])
    chc <- rep(0, length(chu))
    for (i in 1:length(chu)) chc[i] <- sum(z[z[, "chrom"] == 
        chu[i], "weight"])
    for (ich in 1:length(chu)) {
        za <- z[z[, "chrom"] == chu[rev(order(chc))][ich], c("start", 
            "end", "weight"), drop = F]
        za <- za[order(za[, "start"]), , drop = F]
        spairs <- makepairs(za)
        spairs <- spairs[order(spairs[, "start"]), , drop = F]
        psched <- interval.schedule(spairs)
        invl <- 1/(za[, "end"] - za[, "start"] + 1)
        score <- rep(0, nrow(spairs))
        for (sch in 1:max(psched)) {
            pspairs <- spairs[psched == sch, , drop = F]
            ci <- containment.indicator(pspairs[, "start"], pspairs[, 
                "end"], za[, "start"], za[, "end"])
            cigt <- (ci[, 2] >= ci[, 1])
            if (sum(cigt) > 0) {
                pci <- ci[cigt, , drop = F]
                pinvl <- invl[cigt]
                pwa <- za[cigt, "weight"]
                score[psched == sch] <- score[psched == sch] + 
                  (spairs[psched == sch, "end"] - spairs[psched == 
                    sch, "start"] + 1)^pow * fill.values(pci, 
                    pwa * pinvl^pow, nrow(pspairs))
            }
        }
        winloci <- matrix(ncol = 4, nrow = maxmark, dimnames = list(NULL, 
            c("chrom", "start", "end", "score")))
        winloci[, "chrom"] <- chu[rev(order(chc))][ich]
        for (i in 1:maxmark) {
            if (sum(score >= minscore) == 0) 
                break
            spairs <- spairs[score >= minscore, , drop = F]
            score <- score[score >= minscore]
            psched <- interval.schedule(spairs)
            iwin <- which.max(score)
            winloci[i, c("start", "end", "score")] <- c(spairs[iwin, 
                "start"], spairs[iwin, "end"], score[iwin])
            if (ich > 1) 
                if (((sum(accu[, "score"] > winloci[i, 3]) + 
                  i) >= maxmark)) 
                  break
            fixus <- za[, "start"] <= spairs[iwin, "start"] & 
                za[, "end"] >= spairs[iwin, "end"]
            if (sum(fixus) > 0) {
                zaf <- za[fixus, , drop = F]
                finvl <- invl[fixus]
                neweight <- zaf[, "weight"] * (1 - ((spairs[iwin, 
                  "end"] - spairs[iwin, "start"] + 1) * finvl)^pow)
                zaf <- cbind(zaf, neweight)
                for (sch in 1:max(psched)) {
                  pspairs <- spairs[psched == sch, , drop = F]
                  ci <- containment.indicator(pspairs[, "start"], 
                    pspairs[, "end"], zaf[, "start"], zaf[, "end"])
                  cigt <- (ci[, 2] >= ci[, 1])
                  if (sum(cigt) > 0) {
                    pci <- ci[cigt, , drop = F]
                    pinvl <- finvl[cigt]
                    zas <- zaf[cigt, , drop = F]
                    score[psched == sch] <- score[psched == sch] + 
                      (spairs[psched == sch, "end"] - spairs[psched == 
                        sch, "start"] + 1)^pow * fill.values(pci, 
                        pinvl^pow * (zas[, "neweight"] - zas[, 
                          "weight"]), nrow(pspairs))
                  }
                }
                za[fixus, "weight"] <- zaf[, "neweight"]
            }
        }
        winloci <- winloci[!is.na(winloci[, "score"]), , drop = F]
        if (ich == 1) 
            accu <- winloci
        if (ich > 1) {
            accu <- rbind(accu, winloci)
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
