\name{PCORErandomized}
\alias{PCORErandomized}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PCORErandomized(procno = 1, COREobj, boundaries, nprocs = 1, rngoffset = 0)
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
        weight[input[, "chrom"] == coreTable[i, "chrom"]][cigt] <- 0
        weightList[[i]] <- matrix(ncol = 2, data = c(which(input[, 
            "chrom"] == coreTable[i, "chrom"])[cigt], weight[input[, 
            "chrom"] == coreTable[i, "chrom"]][cigt]))
    }
    simscores <- matrix(ncol = nshuffle, nrow = nrow(coreTable))
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
            zaf <- za[za[, "weight"] > 0, , drop = F]
            y <- cbind(c(zaf[, "end"] + 1, zaf[, "start"]), c(-zaf[, 
                "weight"], zaf[, "weight"]))
            y <- y[order(y[, 1]), , drop = F]
            cy2 <- cumsum(y[, 2])
            mscore <- max(cy2)
            for (i in 1:nrow(coreTable)) {
                if (mscore < minscore) {
                  scoremat[i:nrow(coreTable), ich] <- 0
                  break
                }
                scoremat[i, ich] <- mscore
                if (i == nrow(coreTable)) 
                  break
                whereIwas <- which(wza \%in\% weightList[[i]][, 
                  1])
                if (length(whereIwas) > 0) {
                  za[whereIwas, "weight"] <- 0
                  zaf <- za[za[, "weight"] > 0, , drop = F]
                  y <- cbind(c(zaf[, "end"] + 1, zaf[, "start"]), 
                    c(-zaf[, "weight"], zaf[, "weight"]))
                  y <- y[order(y[, 1]), , drop = F]
                  cy2 <- cumsum(y[, 2])
                  mscore <- max(cy2)
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
