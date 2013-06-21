\name{ICORErandomized}
\alias{ICORErandomized}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ICORErandomized(procno = 1, COREobj, boundaries, nprocs = 1, rngoffset = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{procno}{
%%     ~~Describe \code{procno} here~~
}
  \item{COREobj}{
%%     ~~Describe \code{COREobj} here~~
}
  \item{boundaries}{
%%     ~~Describe \code{boundaries} here~~
}
  \item{nprocs}{
%%     ~~Describe \code{nprocs} here~~
}
  \item{rngoffset}{
%%     ~~Describe \code{rngoffset} here~~
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
function (procno = 1, COREobj, boundaries, nprocs = 1, rngoffset = 0) 
{
    need2know <- c("seedme", "nshuffle", "pow", "shufflemethod", 
        "input", "coreTable", "minscore")
    for (item in need2know) assign(item, COREobj[[item]])
    minscore <- max(minscore, min(coreTable[, "score"]))
    set.seed(seedme)
    nshuffle <- nshuffle - rngoffset
    myshuffles <- nshuffle\%/\%nprocs
    shuffleres <- nshuffle\%\%nprocs
    shuffleskip <- myshuffles * (procno - 1) + min(shuffleres, 
        procno - 1)
    if (procno <= shuffleres) 
        myshuffles <- myshuffles + 1
    weightList <- vector(mode = "list", length = nrow(coreTable))
    weight <- input[, "weight"]
    for (i in 1:nrow(coreTable)) {
        za <- input[input[, "chrom"] == coreTable[i, "chrom"], 
            ]
        cigt <- za[, "start"] <= coreTable[i, "start"] & za[, 
            "end"] >= coreTable[i, "end"]
        weight[input[, "chrom"] == coreTable[i, "chrom"]][cigt] <- weight[input[, 
            "chrom"] == coreTable[i, "chrom"]][cigt] * (1 - ((coreTable[i, 
            "end"] - coreTable[i, "start"] + 1)/(za[cigt, "end"] - 
            za[cigt, "start"] + 1))^pow)
        weightList[[i]] <- matrix(ncol = 2, data = c(which(input[, 
            "chrom"] == coreTable[i, "chrom"])[cigt], weight[input[, 
            "chrom"] == coreTable[i, "chrom"]][cigt]))
    }
    simscores <- matrix(ncol = myshuffles, nrow = nrow(coreTable))
    chrmax <- nrow(boundaries)
    advanceRNG(randopt = shufflemethod, nrand = shuffleskip + 
        rngoffset, nevents = nrow(input))
    for (shuffle in 1:myshuffles) {
        if (shufflemethod == "SIMPLE") 
            z <- cbind(randomEventMoves(input[, "end"] - input[, 
                "start"] + 1, boundaries), input[, "weight"])
        if (shufflemethod == "RESCALE") 
            z <- cbind(randomRescaledEventMoves(input[, c("start", 
                "end", "chrom", "weight")], boundaries), input[, 
                "weight"])
        dimnames(z)[[2]] <- c("start", "end", "chrom", "weight")
        chu <- unique(z[, "chrom"])
        chc <- rep(0, length(chu))
        for (i in 1:length(chu)) chc[i] <- sum(z[z[, "chrom"] == 
            chu[i], "weight"])
        scoremat <- matrix(nrow = nrow(coreTable), ncol = length(chu), 
            data = 0)
        for (ich in 1:length(chu)) {
            wza <- which(z[, "chrom"] == chu[rev(order(chc))][ich])
            za <- z[wza, c("start", "end", "weight"), drop = F]
            oza <- order(za[, "start"])
            ooza <- order(oza)
            za <- za[oza, , drop = F]
            spairs <- makepairs(za)
            spairs <- spairs[order(spairs[, "start"]), , drop = F]
            psched <- interval.schedule(spairs)
            invl <- 1/(za[, "end"] - za[, "start"] + 1)
            score <- rep(0, nrow(spairs))
            for (sch in 1:max(psched)) {
                pspairs <- spairs[psched == sch, , drop = F]
                ci <- containment.indicator(pspairs[, "start"], 
                  pspairs[, "end"], za[, "start"], za[, "end"])
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
            for (i in 1:nrow(coreTable)) {
                spairs <- spairs[score >= minscore, , drop = F]
                if (nrow(spairs) == 0) {
                  scoremat[i:nrow(coreTable), ich] <- 0
                  break
                }
                scoremat[i, ich] <- max(score)
                if (i == nrow(coreTable)) 
                  break
                score <- score[score >= minscore]
                psched <- interval.schedule(spairs)
                whereIwas <- which(wza \%in\% weightList[[i]][, 
                  1])
                if (length(whereIwas) > 0) {
                  whereIamNow <- ooza[whereIwas]
                  zaf <- cbind(za[whereIamNow, , drop = F], weightList[[i]][weightList[[i]][, 
                    1] \%in\% wza, 2])
                  dimnames(zaf)[[2]][ncol(zaf)] <- "neweight"
                  za[whereIamNow, "weight"] <- zaf[, "neweight"]
                  zaf <- zaf[order(zaf[, "start"]), , drop = F]
                  invl <- 1/(zaf[, "end"] - zaf[, "start"] + 
                    1)
                  for (sch in 1:max(psched)) {
                    pspairs <- spairs[psched == sch, , drop = F]
                    ci <- containment.indicator(pspairs[, "start"], 
                      pspairs[, "end"], zaf[, "start"], zaf[, 
                        "end"])
                    cigt <- (ci[, 2] >= ci[, 1])
                    if (sum(cigt) > 0) {
                      pci <- ci[cigt, , drop = F]
                      pinvl <- invl[cigt]
                      pwa <- zaf[cigt, "neweight"] - zaf[cigt, 
                        "weight"]
                      score[psched == sch] <- score[psched == 
                        sch] + (spairs[psched == sch, "end"] - 
                        spairs[psched == sch, "start"] + 1)^pow * 
                        fill.values(pci, pwa * pinvl^pow, nrow(pspairs))
                    }
                  }
                }
            }
        }
        simscores[, shuffle] <- apply(scoremat, 1, max)
    }
    return(simscores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
